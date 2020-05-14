# coding: utf-8


# Defines firetasks for parsing Gaussian output files.
import os
import logging
import json
import datetime

from bson.objectid import ObjectId

from pymatgen.io.gaussian import GaussianInput
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from infrastructure.gaussian.utils.utils import get_db, process_run, \
    process_mol, pass_gout_dict, get_chem_schema

logger = logging.getLogger(__name__)

DEFAULT_KEY = 'gout_key'


# TODO: create a util function to do the common things in the PropertiestoDB
#  Firetasks


@explicit_serialize
class ProcessRun(FiretaskBase):
    required_params = ["run"]
    optional_params = ["operation_type", "db", "save_to_db",
                       "save_to_file", "filename", "input_file", "gout_key"]

    def run_task(self, fw_spec):
        run = self["run"]
        operation_type = self.get("operation_type", "get_from_gout")
        input_file = self.get("input_file")
        working_dir = os.getcwd()
        db = self.get('db')

        gout_dict = process_run(operation_type=operation_type, run=run,
                                input_file=input_file, working_dir=working_dir,
                                db=db)
        if "_id" in gout_dict:
            gout_dict["_id"] = ObjectId(gout_dict["_id"])
        if "tag" in fw_spec:
            gout_dict["tag"] = fw_spec["tag"]

        run_list = {}
        if self.get('save_to_db'):
            runs_db = get_db(db)
            run_id = runs_db.insert_run(gout_dict)
            run_list['run_id_list'] = run_id
            logger.info("Saved parsed output to db")

        if self.get('save_to_file'):
            filename = self.get("filename", 'run')
            file = os.path.join(working_dir, f"{filename}.json")
            with open(file, "w") as f:
                f.write(json.dumps(gout_dict, default=DATETIME_HANDLER))
            run_list['run_loc_list'] = file
            logger.info("Saved parsed output to json file")

        uid = self.get("gout_key")
        set_dict = {f"gaussian_output->{DEFAULT_KEY}": gout_dict}
        if uid:
            set_dict[f"gaussian_output->{uid}"] = gout_dict
        # fw_spec = {'gaussian_output: {DEFAULT_KEY: gout_dict, uid: gout_dict}}
        mod_dict = {'_set': set_dict}
        if run_list:
            mod_dict.update({'_push': run_list})
        return FWAction(mod_spec=mod_dict)


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
    optional_params = ["db", "save_to_db", "save_to_file",
                       "additional_prop_doc_fields"]

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
            str(molecules[0].species[index[0]]) + str(index[0]),
            str(molecules[1].species[index[1]]) + str(len(molecules[0]) +
                                                      index[1]))
        be_value = (final_energies[2] -
                    (final_energies[0] + final_energies[1])
                    ) * 27.2114
        mol_schema = get_chem_schema(molecules[2])

        be_dict = {"molecule": molecules[2].as_dict(),
                   "smiles": mol_schema["smiles"],
                   "formula_pretty": mol_schema["formula_pretty"],
                   "energy": final_energies[2],
                   be_key: be_value,
                   "functional": gout_dict[-1]["functional"],
                   "basis": gout_dict[-1]["basis"],
                   "charge": gout_dict[-1]["input"]["charge"],
                   "spin_multiplicity":
                       gout_dict[-1]["input"]["spin_multiplicity"],
                   "phase": gout_dict[-1]["phase"],
                   "tag": gout_dict[-1]["tag"],
                   "state": "successful",
                   "last_updated": datetime.datetime.utcnow()}

        if self.get("additional_prop_doc_fields"):
            be_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            be_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get('save_to_db'):
            db = get_db(db)
            db.insert_property("binding_energy", be_dict,
                               molecules[2].composition.reduced_formula)

        if fw_spec.get("run_loc_list"):
            be_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get('save_to_file'):
            working_dir = os.getcwd()
            if "run_ids" in be_dict:
                del be_dict["run_ids"]
            be_file = os.path.join(working_dir, 'binding_energy.json')
            with open(be_file, "w") as f:
                f.write(json.dumps(be_dict, default=DATETIME_HANDLER))
        logger.info("binding energy calculation complete")
        return FWAction()


@explicit_serialize
class IPEAtoDB(FiretaskBase):
    required_params = ["num_electrons"]
    optional_params = ["db", "save_to_db", "save_to_file",
                       "additional_prop_doc_fields"]

    def run_task(self, fw_spec):
        db = self.get("db")
        num_electrons = self["num_electrons"]
        keys = ["neutral_gas", "anion_gas", "cation_gas",
                "neutral_sol", "anion_sol", "cation_sol"]
        gout_dict = [pass_gout_dict(fw_spec, i) for i in keys]
        molecule = process_mol("get_from_run_dict", gout_dict[0])
        free_energies = [gout['output']['output']["corrections"]
                         ["Gibbs Free Energy"] for gout in gout_dict]
        delta_gibss_sol_neutral = free_energies[3] - free_energies[0]
        delta_gibss_sol_anion = free_energies[4] - free_energies[1]
        delta_gibss_sol_cation = free_energies[5] - free_energies[2]
        delta_gibbs_ox_gas = free_energies[2] - free_energies[0]
        delta_gibbs_red_gas = free_energies[1] - free_energies[0]
        delta_gibbs_ox_sol = delta_gibbs_ox_gas + delta_gibss_sol_cation - \
                             delta_gibss_sol_neutral
        delta_gibbs_red_sol = delta_gibbs_red_gas + delta_gibss_sol_anion - \
                              delta_gibss_sol_neutral
        ea = -delta_gibbs_red_sol * 4.36 * 10 ** -18 * 6.02 * 10 ** 23 / \
             (num_electrons * 9.65 * 10 ** 4)
        ip = -delta_gibbs_ox_sol * 4.36 * 10 ** -18 * 6.02 * 10 ** 23 / \
             (num_electrons * 9.65 * 10 ** 4)

        mol_schema = get_chem_schema(molecule)
        ipea_dict = {"molecule": molecule.as_dict(),
                     "smiles": mol_schema["smiles"],
                     "formula_pretty": mol_schema["formula_pretty"],
                     "IP_ev": ip,
                     "EA_ev": ea,
                     "num_electrons": num_electrons,
                     "functional": gout_dict[0]["functional"],
                     "basis": gout_dict[0]["basis"],
                     "charge": gout_dict[0]["input"]["charge"],
                     "spin_multiplicity":
                         gout_dict[0]["input"]["spin_multiplicity"],
                     "tag": gout_dict[0]["tag"],
                     "state": "successful",
                     "last_updated": datetime.datetime.utcnow()}

        if self.get("additional_prop_doc_fields"):
            ipea_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            ipea_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get('save_to_db'):
            db = get_db(db)
            db.insert_property("ip_ea", ipea_dict,
                               molecule.composition.reduced_formula)

        if fw_spec.get("run_loc_list"):
            ipea_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get('save_to_file'):
            working_dir = os.getcwd()
            if "run_ids" in ipea_dict:
                del ipea_dict["run_ids"]
            be_file = os.path.join(working_dir, 'binding_energy.json')
            with open(be_file, "w") as f:
                f.write(json.dumps(ipea_dict, default=DATETIME_HANDLER))

        logger.info("ip/ea calculation complete")
        return FWAction()
