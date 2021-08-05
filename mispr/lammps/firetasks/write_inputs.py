# coding: utf-8

# Defines firetasks for writing LAMMPS input files

import os
import json
import numpy as np

from string import Template
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.core.structure import Molecule
from pymatgen.io.lammps.data import LammpsData, LammpsDataWrapper
from pymatgen.io.lammps.inputs import write_lammps_inputs
from infrastructure.gaussian.utils.utils import get_mol_formula,\
    get_chem_schema
from infrastructure.lammps.utils.utils import add_ff_labels_to_dict, get_db, \
    process_run
from infrastructure.lammps.database import LammpsSysDb

__author__ = 'Matthew Bliss'
__version__ = '0.0.1'
__email__ = 'matthew.bliss@tufts.edu'
__date__ = 'Apr 14, 2020'

FF_DICT_KEYS = ["Molecule", "Labels", "Masses", "Nonbond", "Bonds", "Angles",
                "Dihedrals", "Impropers", "Improper Topologies", "Charges"]

TEMPLATE_TYPES = ["emin_general", "emin_gaff", "npt", "nvt"]
# ["emin", "npt", "melt", "quench", "nvt"]
# TODO: decide if I want specific templates for melt and quench

TEMPLATE_DIR = os.path.normpath(os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "..", "templates"))

TLEAP_SETTING_DEFAULTS = {"source_file_path": "leaprc.gaff",
                          "mol2_file_path": "mol.mol2",
                          "frcmod_file_path": "mol.frcmod",
                          "prmtop_file_path": "mol.prmtop",
                          "inpcrd_file_path": "mol.inpcrd"}

@explicit_serialize
class WriteDataFile(FiretaskBase):
    _fw_name = "Write Data File"
    required_params = []
    optional_params = ["working_dir", "data_filename", "lammps_data",
                       "lammps_data_wrapper", "system_force_field_dict",
                       "system_mixture_data", "system_box_data",
                       "system_box_data_type", "position_seed",
                       "system_mixture_data_type"]

    def run_task(self, fw_spec):
        # There are three main ways of creating a data file:
        # 1: provide a LammpsData object as "lammps_data" optional parameter
        # 2: provide a LammpsDataWrapper object as "lammps_data_wrapper" \
        #       optional parameter
        # 3: provide the inputs for instantiating a LammpsDataWrapper object

        # Set the path for writing the data file
        working_dir = fw_spec.get("working_dir", self.get("working_dir",
                                                          os.getcwd()))
        os.chdir(working_dir)

        data_file_name = self.get("data_filename", "complex.data")
        data_file_path = os.path.join(working_dir, data_file_name)

        forcefields = None
        mixture = None
        box_data = None
        Lammps_Data_Wrapper = None
        Lammps_Data = None

        # Now to find/instantiate the LammpsData object. There are several \
        #       different input cases
        # Case where LammpsData object is an optional parameter
        if isinstance(self.get("lammps_data"), LammpsData):
            Lammps_Data = self.get("lammps_data")

        # Case where LammpsData object is passed through fw_spec
        elif isinstance(fw_spec.get("lammps_data"), LammpsData):
            Lammps_Data = fw_spec.get("lammps_data")

        # Case where LammpsDataWrapper object is an optional parameter
        elif isinstance(self.get("lammps_data_wrapper"), LammpsDataWrapper):
            Lammps_Data_Wrapper = self.get("lammps_data_wrapper")

        # Case where LammpsDataWrapper object is passed through fw_spec
        elif isinstance(fw_spec.get("lammps_data_wrapper"), LammpsDataWrapper):
            Lammps_Data_Wrapper = fw_spec.get("lammps_data_wrapper")

        # Case where only required arguments for LammpsDataWrapper are provided
        elif all((isinstance(fw_spec.get("system_force_field_dict",
                                         self.get("system_force_field_dict")),
                             dict),
                  isinstance(fw_spec.get("system_mixture_data",
                                         self.get("system_mixture_data")),
                             dict),
                  isinstance(fw_spec.get("system_box_data",
                                         self.get("system_box_data")),
                             (float, int, list, np.ndarray)))):
            forcefields = fw_spec.get("system_force_field_dict",
                                      self.get("system_force_field_dict"))
            # print(forcefields)
            mixture = fw_spec.get("system_mixture_data",
                                  self.get("system_mixture_data"))
            box_data = fw_spec.get("system_box_data",
                                   self.get("system_box_data"))

        # Case where no proper inputs exist: raise error
        else:
            raise KeyError(
                "Not enough input objects for data file creation. "
                "Add optional parameters or check fw_spec."
            )

        # Create LammpsDataWrapper object
        if all((forcefields, mixture, box_data)):
            box_data_type = self.get("system_box_data_type", 'cubic')
            mixture_type = self.get("system_mixture_data_type",
                                    "concentration")
            seed = self.get("position_seed", 150)
            Lammps_Data_Wrapper = LammpsDataWrapper(
                                            system_force_fields=forcefields,
                                            system_mixture_data=mixture,
                                            box_data=box_data,
                                            box_data_type=box_data_type,
                                            mixture_data_type=mixture_type,
                                            seed=seed,
                                            check_ff_duplicates=False)

        # Convert LammpsDataWrapper to LammpsData
        if Lammps_Data_Wrapper:
            Lammps_Data = Lammps_Data_Wrapper.MakeLammpsData()

            n_mols_dict = Lammps_Data_Wrapper.nmol_dict
            smiles_list = [get_chem_schema(forcefields[mol_label]["Molecule"])
                           ["smiles"] for mol_label
                           in Lammps_Data_Wrapper.SortedNames]
            num_mols_list = [n_mols_dict[name] for name in
                             Lammps_Data_Wrapper.SortedNames]
            lmp_box = Lammps_Data.box
            num_atoms_per_mol_list = [
                    len(Lammps_Data_Wrapper.ff_list[name]["Molecule"].sites)
                    for name in Lammps_Data_Wrapper.SortedNames]
            default_masses_list = []
            recalc_masses_list = []
            for name in Lammps_Data_Wrapper.SortedNames:
                default_masses_list += list(Lammps_Data_Wrapper.ff_list[name]
                                            ["Masses"].values())
                for label in Lammps_Data_Wrapper.ff_list[name]["Labels"]:
                    recalc_masses_list.append(
                        Lammps_Data_Wrapper.ff_list[name]["Masses"][label]
                    )



        # Write data file
        if Lammps_Data:
            Lammps_Data.write_file(filename=data_file_path, distance=10,
                                   charge=20)
            return FWAction(update_spec={"smiles": smiles_list,
                                         "nmols": n_mols_dict,
                                         "num_mols_list": num_mols_list,
                                         "num_atoms_per_mol":
                                             num_atoms_per_mol_list,
                                         "default_masses": default_masses_list,
                                         "recalc_masses": recalc_masses_list,
                                         "box": lmp_box})
        else:
            pass

@explicit_serialize
class WriteControlFile(FiretaskBase):
    _fw_name = "Write Control File"
    required_params = []
    optional_params = ["working_dir", "db", "control_filename",
                       "template_filename", "template_dir", "template_str",
                       "control_settings", "save_to_db", "save_to_file"]

    def run_task(self, fw_spec):

        # Set directory for writing control file
        working_dir = fw_spec.get("working_dir", self.get("working_dir",
                                                          os.getcwd()))
        os.makedirs(working_dir, exist_ok=True)
        os.chdir(working_dir)

        db = self.get("db", None)
        save_to_db = self.get("save_to_db", True)
        save_to_file = self.get("save_to_file", False)

        # Set filename for control file
        control_filename = self.get("control_filename", "complex.lammpsin")

        # Set the settings parameter
        control_settings = self.get("control_settings", None)

        # There are three different cases for input of the template:
        # Case 1: template as string
        if isinstance(self.get("template_str"),str):
            template_string = self.get("template_str")
            template_filename = None
        # Cases 2 and 3: template from file, with filename set by \
        #       template_type variable
        elif isinstance(self.get("template_filename"),str):
            template_filename = self.get("template_filename")
            # Case 2: template from file in /../templates/ directory
            if template_filename in TEMPLATE_TYPES:
                template_dir = TEMPLATE_DIR
            # Case 3: template from file provided by user
            elif isinstance(self.get("template_dir"),str):
                template_dir = self.get("template_dir")
            else:
                raise KeyError(
                    "Directory containing custom template file was not " \
                    "specified; add as optional parameter"
                )
        else:
            raise KeyError(
                "Either template was not provided as a valid string or the " \
                "path to a template text file was not provided."
            )

        if template_filename:
            template_path = os.path.join(template_dir, template_filename)
            with open(template_path) as file:
                template_string = file.read()

        write_lammps_inputs(working_dir, template_string,
                            settings=control_settings,
                            script_filename=control_filename)

        smiles_list = fw_spec.get("smiles", [])
        n_mols_dict = fw_spec.get("nmols", {})
        num_mols_list = fw_spec.get("num_mols_list", [])
        lmp_box = fw_spec.get("box", None)
        num_atoms_per_mol_list = fw_spec.get("num_atoms_per_mol", [])
        default_masses_list = fw_spec.get("default_masses", [])
        recalc_masses_list = fw_spec.get("recalc_masses", [])

        run_doc = process_run(smiles_list, n_mols_dict, lmp_box,
                              template_filename, control_settings)
        run_doc.update({"working_dir": working_dir})

        # TODO: add option to not specify db input and still save to database
        #  provided save_to_db is True.
        if db and save_to_db:
            run_db = get_db(db)
            run_db.insert_run(run_doc)
        if save_to_file:
            with open("run_file.json", "w") as run_file:
                json.dump(run_doc, run_file)

        return FWAction(update_spec={"smiles": smiles_list,
                                     "nmols": n_mols_dict,
                                     "num_mols_list": num_mols_list,
                                     "box": lmp_box,
                                     "num_atoms_per_mol":
                                         num_atoms_per_mol_list,
                                     "default_masses": default_masses_list,
                                     "recalc_masses": recalc_masses_list})


# TODO: Write a firetask for writing a general data file
# TODO: Write a firetask for writing a general control file template \
#           (multiple fix options, general template)


@explicit_serialize
class WriteTleapScript(FiretaskBase):
    _fw_name = "Write Tleap Script"
    required_params = []
    optional_params = ["working_dir", "script_string", "template_filename",
                       "template_dir", "tleap_settings", "script_filename"]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        if isinstance(self.get("script_string"), str):
            script_string = self.get("script_string")
        else:
            template_filename = self.get("template_filename", "gaff_tleap")
            template_dir = self.get("template_dir", TEMPLATE_DIR)
            template_file_path = os.path.join(template_dir, template_filename)
            with open(template_file_path) as t_file:
                script_string = t_file.read()

        tleap_settings = TLEAP_SETTING_DEFAULTS.update(
                self.get("tleap_settings", fw_spec.get("tleap_settings", {})))
        # tleap_settings.update(self.get("tleap_settings", \
        #                                   TLEAP_SETTING_DEFAULTS))

        script_string = Template(script_string)\
                                .substitute(tleap_settings,
                                            **TLEAP_SETTING_DEFAULTS)

        script_filename = self.get("script_filename", "tleap.in")
        script_file_path = os.path.join(working_dir, script_filename)

        with open(script_file_path, "w") as script:
            script.write(script_string)


@explicit_serialize
class LabelFFDict(FiretaskBase):
    _fw_name = "Label FF Dict"
    required_params = []
    optional_params = ["mol", "unlabeled_dict", "ff_file", "working_dir",
                       "label", "system_force_field_dict"]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        ff_file = self.get("ff_file")

        if ff_file is not None:
            with open(ff_file, "r") as file:
                unlabeled_dict = json.load(file)
                mol = unlabeled_dict["Molecule"]
        else:

            if not isinstance(self.get("mol"), Molecule):
                raise TypeError('"mol" input must be a PyMatGen Molecule object')

            if not isinstance(self.get("unlabeled_dict"), dict):
                raise TypeError('"unlabeled_dict" input must be a dict')

            mol = self.get("mol")
            unlabeled_dict = self.get("unlabeled_dict")

        if isinstance(self.get("label"), str):
            label = self.get("label", get_mol_formula(mol))

        labeled_dict = add_ff_labels_to_dict(unlabeled_dict, label)

        sys_ff_dict = fw_spec.get("system_force_field_dict",
                                  self.get("system_force_field_dict", {}))
        sys_ff_dict[label] = labeled_dict
        # return FWAction(mod_spec=[{'_set': \
    #                           {"system_force_field_dict": sys_ff_dict}}])
        return FWAction(update_spec={"system_force_field_dict": sys_ff_dict})


@explicit_serialize
class LabelFFDictFromDB(FiretaskBase):
    _fw_name = "Label FF Dict From DB"
    required_params = ["filter"]
    optional_params = ["working_dir", "db", "label", "system_force_field_dict"]

    def run_task(self, fw_spec):
        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        db = self.get("db", None)
        ff_db = get_db(input_db=db)

        filter = self.get("filter")
        unlabeled_ff_dict = ff_db.retrieve_force_field(**filter)
        mol = Molecule.from_dict(unlabeled_ff_dict)
        unlabeled_ff_dict["Molecule"] = mol

        label = self.get("label", get_mol_formula(mol))

        labeled_ff_dict = add_ff_labels_to_dict(unlabeled_ff_dict, label)

        sys_ff_dict = fw_spec.get("system_force_field_dict",
                                  self.get("system_force_field_dict", {}))

        sys_ff_dict[label] = labeled_ff_dict

        return FWAction(update_spec={"system_force_field_dict": sys_ff_dict})


if __name__ == "__main__":
    test_dict = {}
    if isinstance(test_dict.get("gold"),str):
        print(True)
    else:
        print(False)
    gold = test_dict.get("gold")
    print(gold)

    print(os.path.abspath(__file__))
    print(os.path.dirname(os.path.abspath(__file__)))
    print(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "templates"))
    print(os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "templates")))

    s = Template('${who} likes $what')
    d = {'who': 'tim', 'what': 'kung pao'}
    t = {'who': 'george', 'time': 'four o\'clock'}
    p = s.substitute(d, **t)
    print(p)
    print(type(p))