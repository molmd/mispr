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
from infrastructure.gaussian.utils.utils import get_db, process_run, \
    process_mol, pass_gout_dict

logger = logging.getLogger(__name__)

DEFAULT_KEY = 'gout_key'


@explicit_serialize
class ProcessRun(FiretaskBase):
    required_params = ["run"]
    optional_params = ["operation_type", "db", "working_dir", "save_to_db",
                       "save_to_file", "filename", "input_file", "gout_key"]

    def run_task(self, fw_spec):
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
            filename = self.get("filename", 'run')
            file = os.path.join(working_dir, f"{filename}.json")
            with open(file, "w") as f:
                f.write(json.dumps(gout_dict, default=DATETIME_HANDLER))
            logger.info("Saved parsed output to json file")

        uid = self.get("gout_key")
        set_dict = {f"gaussian_output->{DEFAULT_KEY}": gout_dict}
        if uid:
            set_dict[f"gaussian_output->{uid}"] = gout_dict
        # fw_spec = {'gaussian_output: {DEFAULT_KEY: gout_dict, uid: gout_dict}}
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
        # TODO: check if this works in all cases
        #  (if we are saving to db and removing the class)
        # inputs["molecule"] = process_mol("get_from_run_dict", run)
        gaussin = GaussianInput.from_dict(inputs)
        fw_spec["gaussian_input"] = gaussin


@explicit_serialize
class BindingEnergytoDB(FiretaskBase):
    required_params = ["index"]
    optional_params = ["db", "save_to_db", "save_to_file"]

    def run_task(self, fw_spec):
        db = self.get("db")
        index = self["index"]
        keys = ["mol_1", "mol_2", "mol_linked"]
        gout_dict = [pass_gout_dict(fw_spec, i) for i in keys]
        molecules = [process_mol("get_from_run_dict", gout) for gout in
                     gout_dict]
        final_energies = [gout['output']['output']["final_energy"] for gout in
                          gout_dict]
        be_key = "binding_energy_{}_{}_eV".format(
            molecules[0].species[index[0]],
            molecules[1].species[index[1]])
        be_value = (final_energies[2] -
                                (final_energies[0] + final_energies[1])
                                ) * 27.2114

        be_dict = {"molecule": molecules[2].as_dict(),
                   "formula_pretty": molecules[2].composition.reduced_formula,
                   "energy": final_energies[2],
                   be_key: be_value,
                   "state": "successful"}

        if self.get('save_to_db'):
            db = get_db(db)
            db.insert_property("binding_energy", be_dict,
                               molecules[2].composition.reduced_formula)
        if self.get('save_to_file'):
            with open('binding_energy.json', "w") as f:
                f.write(json.dumps(be_dict, default=DATETIME_HANDLER))
        logger.info("binding energy calculation complete")
        return FWAction()

