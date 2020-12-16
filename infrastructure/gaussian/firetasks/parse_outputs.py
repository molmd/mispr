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
    process_mol, pass_gout_dict, get_chem_schema, add_solvent_to_prop_dict

logger = logging.getLogger(__name__)

DEFAULT_KEY = 'gout_key'
HARTREE_TO_EV = 27.2114
HARTREE_TO_KJ = 2600
FARAD = 96.5


# TODO: create a util function to do the common things in the PropertiestoDB
#  Firetasks; ESP and NMR are very similar


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
class NMRtoDB(FiretaskBase):
    required_params = ["keys"]
    optional_params = ["db", "save_to_db", "save_to_file",
                       "additional_prop_doc_fields",
                       "solvent_gaussian_inputs",
                       "solvent_properties"]

    def run_task(self, fw_spec):
        keys = self["keys"]
        db = self.get("db")
        gout_dict = [pass_gout_dict(fw_spec, i) for i in keys]
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = gout_dict + \
                         [pass_gout_dict(fw_spec, i + "_opt") for i in keys]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
        molecule = process_mol("get_from_run_dict", gout_dict[-1])

        mol_schema = get_chem_schema(molecule)
        phase = gout_dict[-1]["phase"]

        # if one calculation is skipped, wall time is considered zero
        run_time = \
            sum([gout.get("wall_time (s)", 0) for gout in full_gout_dict])

        nmr_dict = {"molecule": molecule.as_dict(),
                    "smiles": mol_schema["smiles"],
                    "formula_pretty": mol_schema["formula_pretty"],
                    "energy": gout_dict[-1]['output']['output']["final_energy"],
                    "tensor": gout_dict[-1]["output"]["output"]["tensor"],
                    "functional": gout_dict[-1]["functional"],
                    "basis": gout_dict[-1]["basis"],
                    "charge": gout_dict[-1]["input"]["charge"],
                    "spin_multiplicity":
                        gout_dict[-1]["input"]["spin_multiplicity"],
                    "phase": phase,
                    "tag": gout_dict[-1]["tag"],
                    "state": "successful",
                    "wall_time (s)": run_time,
                    "last_updated": datetime.datetime.utcnow()}

        if phase == "solution":
            solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
            solvent_properties = self.get("solvent_properties", {})
            nmr_dict = add_solvent_to_prop_dict(nmr_dict,
                                                solvent_gaussian_inputs,
                                                solvent_properties)

        if self.get("additional_prop_doc_fields"):
            nmr_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            nmr_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get('save_to_db'):
            db = get_db(db)
            db.insert_property("nmr", nmr_dict,
                               [('formula_pretty', 1), ('smiles', 1)])

        if fw_spec.get("run_loc_list"):
            nmr_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get('save_to_file'):
            working_dir = os.getcwd()
            if "run_ids" in nmr_dict:
                del nmr_dict["run_ids"]
            be_file = os.path.join(working_dir, 'nmr.json')
            with open(be_file, "w") as f:
                f.write(json.dumps(nmr_dict, default=DATETIME_HANDLER))
        logger.info("nmr calculation complete")
        return FWAction()


@explicit_serialize
class BindingEnergytoDB(FiretaskBase):
    required_params = ["index", "keys"]
    optional_params = ["db", "save_to_db", "save_to_file",
                       "additional_prop_doc_fields",
                       "solvent_gaussian_inputs",
                       "solvent_properties"]

    def run_task(self, fw_spec):
        index = self["index"]
        keys = self["keys"]
        db = self.get("db")
        gout_dict = [pass_gout_dict(fw_spec, i) for i in keys]
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = gout_dict + \
                         [pass_gout_dict(fw_spec, i + "_opt") for i in keys]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
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
        phase = gout_dict[-1]["phase"]
        # if one calculation is skipped, wall time is considered zero
        run_time = \
            sum([gout.get("wall_time (s)", 0) for gout in full_gout_dict])

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

        if phase == "solution":
            solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
            solvent_properties = self.get("solvent_properties", {})
            be_dict = add_solvent_to_prop_dict(be_dict,
                                               solvent_gaussian_inputs,
                                               solvent_properties)

        if self.get("additional_prop_doc_fields"):
            be_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            be_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get('save_to_db'):
            db = get_db(db)
            db.insert_property("binding_energy", be_dict,
                               [('formula_pretty', 1), ('smiles', 1)])

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
    required_params = ["num_electrons", "states", "phases", "steps",
                       "root_node_key", "keys", "branch_cation_on_anion"]
    optional_params = ["db", "save_to_db", "save_to_file",
                       "solvent_gaussian_inputs",
                       "solvent_properties",
                       "electrode_potentials",
                       "additional_prop_doc_fields",
                       "gibbs_elec"]

    def run_task(self, fw_spec):
        db = self.get("db")
        num_electrons = self["num_electrons"]
        states = self["states"]
        phases = self["phases"]
        steps = self["steps"]
        root_node_key = self['root_node_key']
        branch_cation_on_anion = self["branch_cation_on_anion"]
        keys = self["keys"]
        gibbs_elec = self.get("gibbs_elec")

        ref_potentials = {'hydrogen': 4.44,
                          'magnesium': 2.07,
                          'lithium': 1.40}
        electrode_potentials = self.get("electrode_potentials") or {}
        electrode_potentials = {
            k.lower(): v for k, v in electrode_potentials.items()}
        electrode_potentials = {**ref_potentials,
                                **electrode_potentials}
        # TODO: get full_gout_dict list including opt to calculate run time
        gout_dict = {i: pass_gout_dict(fw_spec, i) for i in keys}
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = list(gout_dict.values()) + \
                         [pass_gout_dict(fw_spec, i + "_opt") for i in keys]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
        molecule = process_mol("get_from_run_dict",
                               gout_dict[root_node_key])
        final_energies = {i: j['output']['output']["final_energy"]
                          for i, j in gout_dict.items()}
        run_time = sum([gout["wall_time (s)"] for gout in full_gout_dict])

        ip_ea_leafs = {}
        for state in states:
            ip_ea_leafs[state] = {}
            sign = 1 if state == 'anion' else -1
            for phase in phases:
                ip_ea_leafs[state][phase] = {}
                if not branch_cation_on_anion or sign > 0:
                    ip_ea_leafs[state][phase]['full'] = \
                        [f"{phase.lower()}_0_0",
                         f"{phase.lower()}_{sign * num_electrons}_0"]
                    if steps == "multi":
                        for i in range(num_electrons):
                            ip_ea_leafs[state][phase][
                                f'{sign * i}_{sign * (i + 1)}'] = \
                                [f"{phase.lower()}_{sign * i}_0",
                                 f"{phase.lower()}_{sign * (i + 1)}_0"]
                else:
                    ip_ea_leafs[state][phase]['full'] = \
                        [f"{phase.lower()}_0_0",
                         f"{phase.lower()}_0_{num_electrons}"]
                    if steps == "multi":
                        for i in range(num_electrons):
                            ip_ea_leafs[state][phase][f'{i}_{i + 1}'] = \
                                [f"{phase.lower()}_0_{i}",
                                 f"{phase.lower()}_0_{i + 1}"]
        pcet_energy_leafs = {}
        pcet_pka_leafs = {}
        if branch_cation_on_anion:
            for i in range(num_electrons):
                # i = hydrogen if pcet_energy else electron
                for j in range(num_electrons):
                    # j = electron if pcet_energy else hydrogen
                    pcet_energy_leafs[f'{i}_{j}_{j + 1}'] = \
                        [f"solution_{j}_{i}", f"solution_{j + 1}_{i}"]
                    pcet_pka_leafs[f"{i}_{j}_{j + 1}"] = \
                        [f"solution_{i}_{j}", f"solution_{i}_{j + 1}"]

        ip_ea_results = {}
        for phase in phases:
            ip_ea_results[phase] = {}
            for state in states:
                sign = 1 if state == "cation" else -1
                prefix = 'IP' if state == 'cation' else 'EA'
                for k, v in ip_ea_leafs[state][phase].items():
                    res_key = prefix
                    if k != 'full':
                        res_key += '_' + k
                    ip_ea_results[phase][res_key] = \
                        sign * (final_energies[v[1]] - final_energies[v[0]])

                    loc_num_electrons = 1 if k != 'full' else num_electrons
                    if gibbs_elec:
                        if state == "anion":
                            ip_ea_results[phase][res_key] += \
                                gibbs_elec * loc_num_electrons
                        else:
                            ip_ea_results[phase][res_key] -= \
                                gibbs_elec * loc_num_electrons
                    ip_ea_results[phase][res_key] /= \
                        (1 / HARTREE_TO_KJ) * loc_num_electrons * FARAD

        # TODO: cleanup this part
        ip_gas = ip_ea_results.get('gas', {}).get('IP')
        ip_sol = ip_ea_results.get('solution', {}).get('IP')

        ea_gas = ip_ea_results.get('gas', {}).get('EA')
        ea_sol = ip_ea_results.get('solution', {}).get('EA')
        oxidation_sol = {}
        oxidation_gas = {}
        reduction_sol = {}
        reduction_gas = {}
        for key, value in electrode_potentials.items():
            if ip_sol is not None:
                oxidation_sol[key] = ip_sol - value
            if ip_gas is not None:
                oxidation_gas[key] = ip_gas - value
            if ea_sol is not None:
                reduction_sol[key] = ea_sol - value
            if ea_gas is not None:
                reduction_gas[key] = ea_gas - value

        mol_schema = get_chem_schema(molecule)
        ip_ea_dict = {"molecule": molecule.as_dict(),
                      "smiles": mol_schema["smiles"],
                      "formula_pretty": mol_schema["formula_pretty"],
                      "num_electrons": num_electrons,
                      "functional": gout_dict[root_node_key]["functional"],
                      "basis": gout_dict[root_node_key]["basis"],
                      "charge": gout_dict[root_node_key]["input"]["charge"],
                      "spin_multiplicity":
                          gout_dict[root_node_key]["input"][
                              "spin_multiplicity"],
                      "phase": phases,
                      "steps": steps,
                      "gas": {**ip_ea_results['gas'],
                          "electrode_potentials": {
                              "oxidation": oxidation_gas,
                              "reduction": reduction_gas}},
                      "solution": {**ip_ea_results['solution'],
                          "electrode_potentials": {
                              "oxidation": oxidation_sol,
                              "reduction": reduction_sol}},
                      "tag": gout_dict[root_node_key]["tag"],
                      "state": "successful",
                      "wall_time (s)": run_time,
                      "last_updated": datetime.datetime.utcnow()}

        if "solution" in phases:
            solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
            solvent_properties = self.get("solvent_properties", {})
            ip_ea_dict = add_solvent_to_prop_dict(ip_ea_dict,
                                                  solvent_gaussian_inputs,
                                                  solvent_properties)

        if self.get("additional_prop_doc_fields"):
            ip_ea_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            ip_ea_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get('save_to_db'):
            db = get_db(db)
            db.insert_property("ip_ea", ip_ea_dict,
                               ['formula_pretty', 'smiles'])

        if fw_spec.get("run_loc_list"):
            ip_ea_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get('save_to_file'):
            working_dir = os.getcwd()
            if "run_ids" in ip_ea_dict:
                del ip_ea_dict["run_ids"]
            ipea_file = os.path.join(working_dir, 'ip_ea.json')
            with open(ipea_file, "w") as f:
                f.write(json.dumps(ip_ea_dict, default=DATETIME_HANDLER))

        logger.info("ip/ea calculation complete")
        return FWAction()


@explicit_serialize
class BDEtoDB(FiretaskBase):
    required_params = ["bonds", "principle_mol_key", "keys"]
    optional_params = ["db", "save_to_db", "save_to_file",
                       "solvent_gaussian_inputs",
                       "solvent_properties",
                       "additional_prop_doc_fields"]

    def run_task(self, fw_spec):
        db = self.get("db")
        bonds = self["bonds"]
        principle_mol_key = self['principle_mol_key']
        keys = [principle_mol_key] + self["keys"]

        gout_dict = {i: pass_gout_dict(fw_spec, i) for i in keys}
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = list(gout_dict.values()) + \
                         [pass_gout_dict(fw_spec, i + "_opt") for i in keys]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
        molecule = process_mol("get_from_run_dict",
                               gout_dict[principle_mol_key])
        phase = gout_dict[principle_mol_key]["phase"]
        final_energies = {i: j['output']['output']["final_energy"]
                          for i, j in gout_dict.items()}
        run_time = sum([gout["wall_time (s)"] for gout in full_gout_dict])

        fragments = {}
        bde_results = {}
        mol_schema = get_chem_schema(molecule)
        bde_dict = {"molecule": molecule.as_dict(),
                    "smiles": mol_schema["smiles"],
                    "formula_pretty": mol_schema["formula_pretty"],
                    "energy": final_energies[0],
                    "functional": gout_dict[principle_mol_key]["functional"],
                    "basis": gout_dict[principle_mol_key]["basis"],
                    "charge": gout_dict[principle_mol_key]["input"]["charge"],
                    "spin_multiplicity":
                        gout_dict[principle_mol_key]["input"]["spin_multiplicity"],
                    "phase": phase,
                    "fragments": fragments,
                    "bde": bde_results,
                    "tag": gout_dict[principle_mol_key]["tag"],
                    "state": "successful",
                    "wall_time (s)": run_time,
                    "last_updated": datetime.datetime.utcnow()
                    }

        if phase == "solution":
            solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
            solvent_properties = self.get("solvent_properties", {})
            bde_dict = add_solvent_to_prop_dict(bde_dict,
                                                solvent_gaussian_inputs,
                                                solvent_properties)

        if self.get("additional_prop_doc_fields"):
            bde_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            bde_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get('save_to_db'):
            db = get_db(db)
            db.insert_property("nmr", bde_dict,
                               [('formula_pretty', 1), ('smiles', 1)])

        if fw_spec.get("run_loc_list"):
            bde_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get('save_to_file'):
            working_dir = os.getcwd()
            if "run_ids" in bde_dict:
                del bde_dict["run_ids"]
            nmr_file = os.path.join(working_dir, 'nmr.json')
            with open(nmr_file, "w") as f:
                f.write(json.dumps(bde_dict, default=DATETIME_HANDLER))
        logger.info("bde calculation complete")
        return FWAction()
