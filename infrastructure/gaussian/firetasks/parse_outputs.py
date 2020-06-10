# coding: utf-8


# Defines firetasks for parsing Gaussian output files.
import os
import logging
import json
import itertools
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
HARTREE_TO_EV = 27.2114


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
        if "run_time" in fw_spec:
            gout_dict["wall_time (s)"] = fw_spec["run_time"]
        if gout_dict["output"]["has_gaussian_completed"]:
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
        else:
            raise ValueError(f"Gaussian did not complete normally, Terminating")


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
            run = pass_gout_dict(fw_spec, DEFAULT_KEY)
            # run = fw_spec["gaussian_output"][DEFAULT_KEY]

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
                    ) * HARTREE_TO_EV
        mol_schema = get_chem_schema(molecules[2])
        # if one calculation is skipped, wall time is considered zero
        run_time = sum([gout.get("wall_time (s)", 0) for gout in gout_dict])

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
                   "wall_time (s)": run_time,
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
    # TODO: cleanup this firetask without assuming we have specific calcs
    required_params = ["num_electrons", "states", "phases", "parents_dict"]
    optional_params = ["db", "save_to_db", "save_to_file",
                       "solvent_gaussian_inputs",
                       "solvent_properties",
                       "electrode_potentials",
                       "additional_prop_doc_fields"]

    def run_task(self, fw_spec):
        db = self.get("db")
        num_electrons = self["num_electrons"]
        states = self["states"]
        phases = self["phases"]
        parents_dict = self["parents_dict"]
        if "solution" in phases:
            if not self.get("solvent_gaussian_inputs"):
                solvent = "water"
            else:
                solvent_inputs = [
                    i.lower() for i in
                    self.get("solvent_gaussian_inputs").strip("()").split(",")]
                solvent = [
                    string for string in solvent_inputs if "solvent" in
                                                           string][0].split("=")[1]
                solvent_properties = self.get("solvent_properties", {})

        ref_potentials = {'hydrogen': 4.44,
                          'magnesium': 2.07,
                          'lithium': 1.40}
        electrode_potentials = {
            k.lower(): v for k, v in
            self.get("electrode_potentials", {}).items()}
        electrode_potentials = {**ref_potentials,
                                **electrode_potentials}

        delta_g = {}
        keys = ["{}_{}".format(i.lower(), j.lower())
                for i, j in itertools.product(states, phases)]
        gout_dict = {i: pass_gout_dict(fw_spec, i) for i in keys}
        molecule = process_mol("get_from_run_dict", gout_dict[0])
        final_energies = {i: j['output']['output']["final_energy"]
                          for i, j in gout_dict}
        run_time = sum([gout["wall_time (s)"] for gout in gout_dict])
        for i, j in parents_dict.items():
            delta_g[i + '-' + j] = final_energies[i] - final_energies[j]
        for phase in phases:
            filtered_states = [i.lower() for i in states if i.lower()
                               in ['anion', 'cation']]
            for state in filtered_states:
                key_1 = "{}_{}".format(state, phase.lower())
                key_2 = "reference_{}".format(phase.lower())
                key = key_1 + '-' + key_2
                if key in delta_g:
                    continue
                delta_g[key] = 0
                for i, sign_i in zip([key_1, key_2], [1, -1]):
                    current = i
                    while current:
                        new = parents_dict.get(current)
                        if new:
                            delta_g[key] += sign_i * delta_g[
                                current + '-' + new]
                        current = new

        # TODO: cleanup this part
        ip_gas = delta_g[
                     "cation_gas-reference_gas"] * HARTREE_TO_EV / num_electrons
        ip_sol = delta_g[
                     "cation_solution-reference_solution"] * HARTREE_TO_EV / num_electrons
        ea_sol = -delta_g[
            "anion_solution-reference_solution"] * HARTREE_TO_EV / num_electrons
        ea_gas = -delta_g[
            "anion_gas-reference_gas"] * HARTREE_TO_EV / num_electrons
        oxidation_sol = {}
        oxidation_gas = {}
        reduction_sol = {}
        reduction_gas = {}
        for key, value in electrode_potentials.items():
            oxidation_sol[key] = ip_sol - value
            oxidation_gas[key] = ip_gas - value
            reduction_sol[key] = ea_sol - value
            reduction_gas[key] = ea_gas - value

        mol_schema = get_chem_schema(molecule)
        ipea_dict = {"molecule": molecule.as_dict(),
                     "smiles": mol_schema["smiles"],
                     "formula_pretty": mol_schema["formula_pretty"],
                     "num_electrons": num_electrons,
                     "functional": gout_dict[0]["functional"],
                     "basis": gout_dict[0]["basis"],
                     "charge": gout_dict[0]["input"]["charge"],
                     "spin_multiplicity":
                         gout_dict[0]["input"]["spin_multiplicity"],
                     "phase": phases,
                     "solvent": solvent,
                     "solvent_properties": solvent_properties,
                     "solution": {"IP": ip_sol,
                                  "EA": ea_sol,
                                  "electrode_potentials": {
                                      "oxidation": oxidation_sol,
                                      "reduction": reduction_sol}
                                  },
                     "gas": {"IP": ip_gas,
                             "EA": ea_gas,
                             "electrode_potentials": {
                                 "oxidation": oxidation_gas,
                                 "reduction": reduction_gas}
                             },
                     "tag": gout_dict[0]["tag"],
                     "state": "successful",
                     "wall_time (s)": run_time,
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
            ipea_file = os.path.join(working_dir, 'ip_ea.json')
            with open(ipea_file, "w") as f:
                f.write(json.dumps(ipea_dict, default=DATETIME_HANDLER))

        logger.info("ip/ea calculation complete")
        return FWAction()
