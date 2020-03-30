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

logger = logging.getLogger(__name__)

DEFAULT_KEY = 'gout_key'


@explicit_serialize
class ProcessRun(FiretaskBase):
    required_params = ["run"]
    optional_params = ["operation_type", "db", "working_dir", "save_to_db",
                       "save_to_file", "filename", "input_file", "gout_key"]

    def run_task(self, fw_spec):
        # TODO: prevent redundant saving (if operation type is from db)
        run = self["run"]
        operation_type = self.get("operation_type", "get_from_gout")
        input_file = self.get("input_file")
        working_dir = self.get('working_dir', os.getcwd())
        db = self.get('db')

        gout_dict = process_run(operation_type=operation_type, run=run,
                                input_file=input_file, working_dir=working_dir,
                                db=db)
        if "_id" in gout_dict:
            gout_dict["_id"] = ObjectId(gout_dict["_id"])
        if "tag" in fw_spec:
            gout_dict["tag"] = fw_spec["tag"]

        if self.get('save_to_db'):
            runs_db = get_db(db)
            runs_db.insert_run(gout_dict)
            logger.info("Saved parsed output to db")
        if self.get('save_to_file'):
            with open(self.get('filename', 'run.json'), "w") as f:
                f.write(json.dumps(gout_dict, default=DATETIME_HANDLER))
            logger.info("Saved parsed output to json file")

        uid = self.get("gout_key")
        set_dict = {f"gaussian_output->{DEFAULT_KEY}": gout_dict}
        if uid:
            set_dict[f"gaussian_output->{uid}"] = gout_dict
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
            inputs[f'{k}'] = self.get("gaussian_input_params", {}). \
                get(f'{k}', run['input'].get(f'{k}'))
        inputs['molecule'] = run['output']['output']['molecule']
        gaussin = GaussianInput.from_dict(inputs)
        fw_spec["gaussian_input"] = gaussin


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