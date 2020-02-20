# coding: utf-8


# Defines firetasks for parsing Gaussian output files.
import os
import logging
import json

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from infrastructure.gaussian.utils.utils import get_chem_schema

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
    """
    required_params = ["db"]
    optional_params = [
        "working_dir", "input_file", "output_file", "fw_spec_field"
    ]

    def run_task(self, fw_spec):
        from pymatgen.io.gaussian import GaussianOutput, GaussianInput
        from pymatgen.core.structure import Molecule
        from infrastructure.gaussian.database import GaussianCalcDb
        # TODO: Include multirun format
        # TODO: cleanup task_doc
        working_dir = self.get("working_dir", os.getcwd())
        input_file = self.get("input_file", "mol.com")
        output_file = self.get("output_file", "mol.out")

        logger.info("Parsing directory: {}".format(working_dir))

        output_path = os.path.join(working_dir, output_file)
        gout = GaussianOutput(output_path).as_dict()
        gout['input']['charge'] = gout['charge']
        gout['input']['spin_multiplicity'] = gout['spin_multiplicity']

        if self.get("input_file"):
            input_path = os.path.join(working_dir, input_file)
            gin = GaussianInput.from_file(input_path).as_dict()

            gout['input']['input_parameters'] = gin['input_parameters']
            task_doc = {'input': gin, 'output': gout}
        else:
            task_doc = {'output': gout}
            logger.info("Input parameters at the end of the Gaussian input "
                        "section will not be saved to the database due to a "
                        "missing input file")

        if "tag" in fw_spec:
            task_doc.update({"tag": fw_spec["tag"]})

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update({self.get("fw_spec_field"):
                                 fw_spec.get(self.get("fw_spec_field"))})

        task_doc.update(
            get_chem_schema(Molecule.from_dict(gout['output']['molecule']))
        )

        task_doc['functional'] = gout['input']['functional']
        task_doc['basis'] = gout['input']['basis_set']
        job_types = \
            list(filter(lambda x: x in
                                  {k.lower(): v for k, v in
                                   gout['input']['route_parameters'].items()},
                        JOB_TYPES))
        task_doc['type'] = ';'.join(job_types)
        task_doc = json.loads(json.dumps(task_doc))
        if isinstance(self.get('db'), dict):
            runs_db = GaussianCalcDb(**self['db'])
        else:
            runs_db = GaussianCalcDb.from_db_file(self.get('db'))
        gout_id = runs_db.insert_run(task_doc)
        logger.info("Finished parsing output and saving to db")
        condition = True
        if not fw_spec.get('additional_fields'):
            condition = False
        else:
            for i in fw_spec['additional_fields']:
                condition = condition and task_doc['output'][i]

        return FWAction(update_spec={"gaussian_output_id": gout_id})
