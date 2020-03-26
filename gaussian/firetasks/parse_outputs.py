# coding: utf-8


# Defines firetasks for parsing Gaussian output files.
import os
import logging
import json

from bson.objectid import ObjectId

from pymatgen.io.gaussian import GaussianInput
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from infrastructure.gaussian.utils.utils import get_chem_schema, get_db, \
    get_run_from_fw_spec, process_run
from infrastructure.gaussian.utils.utils import recursive_signature_remove


logger = logging.getLogger(__name__)

JOB_TYPES = {'sp', 'opt', 'freq', 'irc', 'ircmax', 'scan', 'polar', 'admp',
             'bomd', 'eet', 'force', 'stable', 'volume', 'density', 'guess',
             'pop', 'scrf', 'cphf', 'prop', 'nmr', 'cis', 'zindo', 'td', 'eom',
             'sac-ci'}
DEFAULT_KEY = 'gout_key'


@explicit_serialize
class GaussianToDB(FiretaskBase):
    """
    Assumes the functional/basis_set to be the first parameters on the
    route_parameters line.
    tag is required as an fw_spec
    """
    required_params = []
    optional_params = ["db", "working_dir", "input_file", "output_file",
                       "fw_spec_field", "save_mol_file", "fmt", "filename",
                       "gout_id_key"]

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
        mol = Molecule.from_dict(gout['output']['molecule'])
        task_doc.update(get_chem_schema(mol))
        task_doc['functional'] = task_doc['input']['functional']
        task_doc['basis'] = task_doc['input']['basis_set']
        job_types = \
            list(filter(lambda x: x in
                                  {k.lower(): v for k, v in
                                   task_doc['input'][
                                       'route_parameters'].items()},
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

        uid = self.get("gout_id_key")
        set_dict = {f"gaussian_output_id->{DEFAULT_KEY}": gout_id}
        if uid:
            set_dict[f"gaussian_output_id->{uid}"] = gout_id

        return FWAction(mod_spec={'_set': set_dict})


@explicit_serialize
class BindingEnergytoDB(FiretaskBase):
    required_params = ["keys", 'main_run_key', 'new_prop']
    optional_params = ["db"]

    def run_task(self, fw_spec):
        run_db = get_db(self.get('db'))
        runs = [get_run_from_fw_spec(fw_spec, i, run_db) for i in self['keys']]
        props = [i['output']['output']["final_energy"] for i in runs]
        result = (props[2] - (props[0] + props[1])) * 27.2114
        main_run = get_run_from_fw_spec(fw_spec, self['main_run_key'], run_db)
        run_db.update_run(new_values={self["new_prop"]: result},
                          _id=ObjectId(main_run['_id']))
        logger.info("binding energy calculation complete")


@explicit_serialize
class ProcessRun(FiretaskBase):
    required_params = ["run"]
    optional_params = ["operation_type", "db", "working_dir", "save_to_db",
                       "save_to_file", "filename", "input_file", "gout_key"]

    @staticmethod
    def _modify_gout(gout):
        gout['input']['charge'] = gout['charge']
        gout['input']['spin_multiplicity'] = gout['spin_multiplicity']
        del_keys_out = ('nsites', 'unit_cell_formula',
                        'reduced_cell_formula', 'pretty_formula',
                        'elements', 'nelements', 'charge',
                        'spin_multiplicity')
        [gout.pop(k, None) for k in del_keys_out]

    @staticmethod
    def _create_gin(gout, working_dir, input_file):
        from pymatgen.io.gaussian import GaussianInput
        if input_file:
            input_path = os.path.join(working_dir, input_file)
            gin = GaussianInput.from_file(input_path).as_dict()
            gin['nbasisfunctions'] = gout['input']['nbasisfunctions']
            gin['pcm_parameters'] = gout['input']['pcm_parameters']
            return gin
        else:
            gin = gout['input']
            gin['input_parameters'] = None
            gin['@class'] = 'GaussianInput'
            gin['@module'] = 'pymatgen.io.gaussian'
            logger.info("Input parameters at the end of the Gaussian input "
                        "section will not be saved to the database due to "
                        "a missing input file")
            return gin

    @staticmethod
    def _job_types(gin):
        return list(filter(lambda x: x in {k.lower(): v for k, v in
                                           gin['route_parameters'].items()},
                           JOB_TYPES))

    def run_task(self, fw_spec):
        from pymatgen.core.structure import Molecule

        run = self["run"]
        operation_type = self.get("operation_type", "get_from_gout")
        working_dir = self.get('working_dir', os.getcwd())
        db = self.get('db')

        gout = process_run(operation_type=operation_type, run=run,
                           working_dir=working_dir, db=db)
        task_doc = {'output': gout}
        if self.get('save_to_db') or self.get('save_to_file'):
            # TODO: check if other solvation models are supported
            self._modify_gout(gout)
            gin = self._create_gin(gout, working_dir, self.get("input_file"))
            del gout['input']
            job_types = self._job_types(gin)
            mol = Molecule.from_dict(gout['output']['molecule'])

            task_doc = {'input': gin, 'output': gout, 'tag': fw_spec["tag"],
                        'functional': gin['functional'],
                        'basis': gin['basis_set'],
                        'phase': 'solution' if gout['is_pcm'] else 'gas',
                        'type': ';'.join(job_types),
                        **get_chem_schema(mol)}
            task_doc = {i: j for i, j in task_doc.items() if i not in
                        ['sites', '@module', '@class', 'charge',
                         'spin_multiplicity']}
        task_doc = json.loads(json.dumps(task_doc))

        if self.get('save_to_db'):
            runs_db = get_db(db)
            runs_db.insert_run(task_doc)
            logger.info("Saved parsed file to db")
        if self.get('save_to_file'):
            with open(self.get('filename', 'run.json'), "w") as f:
                f.write(json.dumps(task_doc, default=DATETIME_HANDLER))
            logger.info("Saved parsed output to json file")

        task_doc = recursive_signature_remove(task_doc)

        uid = self.get("gout_key")
        set_dict = {f"gaussian_output->{DEFAULT_KEY}": task_doc}
        if uid:
            set_dict[f"gaussian_output->{uid}"] = task_doc
        return FWAction(mod_spec={'_set': set_dict})


@explicit_serialize
class RetrieveGaussianOutput(FiretaskBase):
    """
    Returns a Gaussian output object from the database and converts it to a
    Gaussian input object
    """
    required_params = []
    optional_params = ["db", "gaussian_input_params", "run_id", "smiles",
                       "functional", "basis", "type", "phase", "tag"]

    def run_task(self, fw_spec):

        # if a Gaussian output dict is passed through fw_spec
        if fw_spec.get("gaussian_output"):
            run = fw_spec["gaussian_output"][DEFAULT_KEY]

        elif self.get("run_id"):
            run = process_run(operation_type="get_from_run_id",
                              run=self.get("run_id"), db=self.get("db"))

        # if a Gaussian output dictionary is retrieved from db
        else:
            query = {'smiles': self.get('smiles'), 'type': self.get('type'),
                     'functional': self.get('functional'),
                     'basis': self.get('basis'), 'phase': self.get('phase')}
            if 'tag' in self:
                query['tag'] = self['tag']
            run = process_run(operation_type="get_from_run_query", run=query,
                              db=self.get('db'))

        # create a gaussian input object from run
        if self.get("gaussian_input_params") is None:
            logger.info("No gaussian input parameters provided; will use "
                        "run parameters")
        inputs = {}
        for k, v in run['input'].items():
            # use gaussian_input_params if defined, otherwise use run parameters
            inputs[f'{k}'] = self.get("gaussian_input_params", {}).\
                get(f'{k}', run['input'].get(f'{k}'))
        inputs['molecule'] = run['output']['output']['molecule']
        gaussin = GaussianInput.from_dict(inputs)
        fw_spec["gaussian_input"] = gaussin
