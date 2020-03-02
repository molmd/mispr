# coding: utf-8


# Defines firetasks for parsing Gaussian output files.
import os
import logging
import json

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from infrastructure.gaussian.utils.utils import get_chem_schema, get_db

logger = logging.getLogger(__name__)


JOB_TYPES = {'sp', 'opt', 'freq', 'irc', 'ircmax', 'scan', 'polar', 'admp',
             'bomd', 'eet', 'force', 'stable', 'volume', 'density', 'guess',
             'pop', 'scrf', 'cphf', 'prop', 'nmr', 'cis', 'zindo', 'td', 'eom',
             'sac-ci'}


@explicit_serialize
class GaussianToDB(FiretaskBase):
    """
    Assumes the functional/basis_set to be the first parameters on the
    route_parameters line.
    tag is required as an fw_spec
    """
    required_params = []
    optional_params = ["db", "working_dir", "input_file", "output_file",
                       "fw_spec_field"]

    def run_task(self, fw_spec):
        from pymatgen.io.gaussian import GaussianOutput, GaussianInput
        from pymatgen.core.structure import Molecule
        from infrastructure.gaussian.database import GaussianCalcDb
        # TODO: Include multirun format
        working_dir = self.get("working_dir", os.getcwd())
        input_file = self.get("input_file", "mol.com")
        output_file = self.get("output_file", "mol.out")

        logger.info("Parsing directory: {}".format(working_dir))

        output_path = os.path.join(working_dir, output_file)
        gout = GaussianOutput(output_path).as_dict()
        gout['input']['charge'] = gout['charge']
        gout['input']['spin_multiplicity'] = gout['spin_multiplicity']
        del_keys_out = ('nsites', 'unit_cell_formula', 'reduced_cell_formula',
                        'pretty_formula', 'elements', 'nelements', 'charge',
                        'spin_multiplicity')
        [gout.pop(k, None) for k in del_keys_out]
        # TODO: check if other solvation models are supported
        if gout['is_pcm']:
            phase = 'solution'
        else:
            phase = 'gas'
        if self.get("input_file"):
            input_path = os.path.join(working_dir, input_file)
            gin = GaussianInput.from_file(input_path).as_dict()
            gin['nbasisfunctions'] = gout['input']['nbasisfunctions']
            gin['pcm_parameters'] = gout['input']['pcm_parameters']
            del gout['input']
            task_doc = {'input': gin, 'output': gout}
        else:
            task_doc = {'output': gout}
            task_doc['input'] = task_doc['output']['input']
            task_doc['input']['input_parameters'] = None
            task_doc['input']['@class'] = 'GaussianInput'
            task_doc['input']['@module'] = 'pymatgen.io.gaussian'
            del task_doc['output']['input']
            logger.info("Input parameters at the end of the Gaussian input "
                        "section will not be saved to the database due to a "
                        "missing input file")

        task_doc.update({"tag": fw_spec["tag"]})
        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update({self.get("fw_spec_field"):
                                 fw_spec.get(self.get("fw_spec_field"))})
        task_doc.update(
            get_chem_schema(Molecule.from_dict(gout['output']['molecule']))
        )
        task_doc['functional'] = task_doc['input']['functional']
        task_doc['basis'] = task_doc['input']['basis_set']
        job_types = \
            list(filter(lambda x: x in
                                  {k.lower(): v for k, v in
                                   task_doc['input']['route_parameters'].items()},
                        JOB_TYPES))
        task_doc['type'] = ';'.join(job_types)
        task_doc['phase'] = phase
        del_keys_doc = ('sites', '@module', '@class', 'charge',
                        'spin_multiplicity')
        [task_doc.pop(k, None) for k in del_keys_doc]
        task_doc = json.loads(json.dumps(task_doc))

        runs_db = get_db(self.get('db'))

        gout_id = runs_db.insert_run(task_doc)
        logger.info("Finished parsing output and saving to db")

        return FWAction(update_spec={"gaussian_output_id": gout_id})
