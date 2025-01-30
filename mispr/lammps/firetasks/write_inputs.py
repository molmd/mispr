"""Define firetasks for writing LAMMPS input files."""

import os
import json
import logging

from string import Template

import numpy as np

from pymatgen.core.structure import Molecule
from pymatgen.io.lammps.data import LammpsData, LammpsDataWrapper
from pymatgen.io.lammps.inputs import write_lammps_inputs

from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from mispr.lammps.defaults import TEMPLATE_TYPES, TLEAP_SETTINGS
from mispr.lammps.utilities.utilities import (
    get_db,
    process_run,
    add_ff_labels_to_dict,
    lammps_mass_to_element,
)
from mispr.gaussian.utilities.metadata import get_chem_schema, get_mol_formula

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.4"

logger = logging.getLogger(__name__)

TEMPLATE_DIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "templates")
)
DEFAULT_KEY = "lammpsin_key"


@explicit_serialize
class WriteDataFile(FiretaskBase):
    """
    Write LAMMPS data file for a bulk system in the current working directory. User can
    provide a ``LammpsData object``, a ``LammpsDataWrapper`` object, or the inputs for
    instantiating a ``LammpsDataWrapper`` object. See the ``LammpsData`` and
    ``LammpsDataWrapper`` classes in the ``pymatgen.io.lammps.data`` module for more
    information.

    Args:
        working_dir (str, optional): Path to the working directory. Defaults to the
            current working directory.
        data_filename (str, optional): Name of data file to be written. Defaults to
            "complex.data".
        lammps_data (LammpsData, optional): The ``LammpsData`` object to be written as
            the data file. Can be passed through the ``fw_spec`` or input by the user.
            If both are provided, the one provided directly will be used.
        lammps_data_wrapper (LammpsDataWrapper, optional): The wrapper for the
            ``LammpsData`` object to be written as data file. Can be passed through
            ``fw_spec`` or input by the user. If both are provided, the one provided
            directly will be used. Will be ignored if ``lammps_data`` is provided.
        system_force_field_dict (dict, optional): The force field parameters for the
            system. Used as input for ``LammpsDataWrapper``. Intended to be created by
            the ``GetForceField`` task. Can be passed through ``fw_spec`` or input by
            the user. If both are provided, will take the one provided in the ``fw_spec``.
        system_mixture_data (dict, optional): The information that will be used to
            calculate the number of molecules of each type in the system. Used as input
            for a ``LammpsDataWrapper`` object. Can be read from ``fw_spec`` or input
            by the user; If both are provided, will use the one provided in the ``fw_spec``.
        system_box_data (float, int, list, ndarray, optional): Information about the
            system box; the value of this arg is determined by the
            ``system_box_data_type``. Used as input for a ``LammpsDataWrapper`` object.
            Can be passed through ``fw_spec`` or input by the user. If both are provided,
            will take the one provided in the ``fw_spec``.
        system_box_data_type (str, optional): Determines the type for ``system_box_data``.
            Used as input for a ``LammpsDataWrapper`` object. Defaults to "cubic".
        position_seed (int, optional): Seed for randomizing the positions of atoms.
            Indirectly used as input for packmol through ``LammpsDataWrapper``. Defaults
            to 150.
        system_mixture_data_type (str, optional): Determines the content of
            ``system_mixture_data``. Used as input for a ``LammpsDataWrapper`` object.
            Defaults to "concentration".
        scale_charges (bool, optional): Whether to scale the partial charges of all
            molecules that have a non-zero net charge in the system. Only used if
            building a system from ``LammpsDataWrapper`` inputs.
        charge_scaling_factor (float): Factor by which to scale charges in the system.
            Only used if building a system from ``LammpsDataWrapper`` inputs and if
            ``scale_charges`` is ``True``.
    """
    _fw_name = "Write Data File"
    required_params = []
    optional_params = [
        "working_dir",
        "data_filename",
        "lammps_data",
        "lammps_data_wrapper",
        "system_force_field_dict",
        "system_mixture_data",
        "system_box_data",
        "system_box_data_type",
        "position_seed",
        "system_mixture_data_type",
        "scale_charges",
        "charge_scaling_factor",
    ]

    def run_task(self, fw_spec):
        # There are three main ways of creating a data file:
        # 1: provide a LammpsData object as "lammps_data"
        #       optional parameter
        # 2: provide a LammpsDataWrapper object as
        #       "lammps_data_wrapper" optional parameter
        # 3: provide the inputs for instantiating a LammpsDataWrapper
        #       object

        # Set the path for writing the data file
        working_dir = fw_spec.get("working_dir",
                                  self.get("working_dir", os.getcwd()))
        os.chdir(working_dir)

        data_file_name = self.get("data_filename", "complex.data")
        data_file_path = os.path.join(working_dir, data_file_name)

        force_fields = None
        mixture = None
        box_data = None
        lammps_data_wrapper = None
        lammps_data = None

        # Now to find/instantiate the LammpsData object. There are
        #   several different input cases:
        # Case where LammpsData object is an optional parameter
        if isinstance(self.get("lammps_data"), LammpsData):
            lammps_data = self.get("lammps_data")

        # Case where LammpsData object is passed through fw_spec
        elif isinstance(fw_spec.get("lammps_data"), LammpsData):
            lammps_data = fw_spec.get("lammps_data")

        # Case where LammpsDataWrapper object is an optional parameter
        elif isinstance(self.get("lammps_data_wrapper"), LammpsDataWrapper):
            lammps_data_wrapper = self.get("lammps_data_wrapper")

        # Case where LammpsDataWrapper object is passed through fw_spec
        elif isinstance(fw_spec.get("lammps_data_wrapper"), LammpsDataWrapper):
            lammps_data_wrapper = fw_spec.get("lammps_data_wrapper")

        # Case where only required arguments for LammpsDataWrapper
        #   are provided
        elif all(
            (
                isinstance(
                    fw_spec.get("system_force_field_dict",
                                self.get("system_force_field_dict")),
                    dict,
                ),
                isinstance(
                    fw_spec.get("system_mixture_data",
                                self.get("system_mixture_data")),
                    dict,
                ),
                isinstance(
                    fw_spec.get("system_box_data",
                                self.get("system_box_data")),
                    (float, int, list, np.ndarray),
                ),
            )
        ):
            force_fields = fw_spec.get("system_force_field_dict",
                                       self.get("system_force_field_dict"))
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
        if all((force_fields, mixture, box_data)):
            scaling_factor = self.get("charge_scaling_factor")
            if self.get("scale_charges") and scaling_factor:
                for k, v in force_fields.items():
                    if round(sum(force_fields[k]["Charges"])) != 0:
                        force_fields[k]["Charges"] = [
                            i * scaling_factor for i in force_fields[k]["Charges"]
                        ]
            box_data_type = self.get("system_box_data_type", "cubic")
            mixture_type = self.get("system_mixture_data_type",
                                    "concentration")
            seed = self.get("position_seed", 150)
            lammps_data_wrapper = LammpsDataWrapper(
                system_force_fields=force_fields,
                system_mixture_data=mixture,
                box_data=box_data,
                box_data_type=box_data_type,
                mixture_data_type=mixture_type,
                seed=seed,
                check_ff_duplicates=False,
            )

        # Convert LammpsDataWrapper to LammpsData
        if lammps_data_wrapper:
            lammps_data = lammps_data_wrapper.build_lammps_data()

            n_mols_dict = lammps_data_wrapper.nmol_dict
            molecules_list = [
                force_fields[mol_label]["Molecule"]
                for mol_label in lammps_data_wrapper.sorted_mol_names
            ]
            smiles_list = [
                get_chem_schema(force_fields[mol_label]["Molecule"])["smiles"]
                for mol_label in lammps_data_wrapper.sorted_mol_names
            ]
            num_mols_list = [
                n_mols_dict[name]
                for name in lammps_data_wrapper.sorted_mol_names
            ]
            lmp_box = lammps_data.box
            num_atoms_per_mol_list = [
                len(lammps_data_wrapper.ff_list[name]["Molecule"].sites)
                for name in lammps_data_wrapper.sorted_mol_names
            ]
            default_masses_list = []
            recalc_masses_list = []
            for name in lammps_data_wrapper.sorted_mol_names:
                default_masses_list += list(
                    lammps_data_wrapper.ff_list[name]["Masses"].values()
                )
                for label in lammps_data_wrapper.ff_list[name]["Labels"]:
                    recalc_masses_list.append(
                        lammps_data_wrapper.ff_list[name]["Masses"][label]
                    )

        # Write data file
        if lammps_data:
            lammps_data.write_file(filename=data_file_path,
                                   distance=10,
                                   charge=20)
            return FWAction(
                update_spec={
                    "smiles": smiles_list,
                    "nmols": n_mols_dict,
                    "num_mols_list": num_mols_list,
                    "num_atoms_per_mol": num_atoms_per_mol_list,
                    "default_masses": default_masses_list,
                    "recalc_masses": recalc_masses_list,
                    "molecules": molecules_list,
                    "box": lmp_box,
                },
                propagate=True,
            )
        else:
            pass


@explicit_serialize
class WriteControlFile(FiretaskBase):
    """
    Write a LAMMPS control file based on a template file or string.

    Args:
        working_dir (str, optional): Path to working directory. Default is current
            working directory.
        db (str or dict, optional): Connection information for the calc database; If
            ``save_runs_to_db`` is ``True``, then a document will be stored in the
            database. Default is ``None``.
        control_filename (str, optional): Name of control file to be written. Default is
            "complex.lammpsin".
        template_filename (str, optional): Name of template control file to be used.
            Is not used if ``template_str`` is provided. Can refer to a custom template
            or a template name from the ``mispr/lammps/templates`` directory.
        template_dir (str, optional): Path to directory containing template file. If the
            template is from ``mispr/lammps/templates``, then this is not needed.
        template_str (str, optional): String containing template for control file.
        control_settings (dict, optional): Dictionary of settings to be used in the
            control file. If a mispr template is used, this dictionary will update the
            corresponding control_settings dictionary in the ``mispr/lammps/defaults.py``
            file. Default is None.
        save_runs_to_db (bool, optional): Whether to save the run info to the calc
            database. Default is ``False``.
        save_runs_to_file (bool, optional): Whether to save the run info to a file in
            the working directory. Default is ``True``.
        lammpsin_key (str, optional): Key in the ``fw_spec`` that corresponds with the
            ``run_doc`` information. Default is "lammpsin".
    """
    _fw_name = "Write Control File"
    required_params = []
    optional_params = [
        "working_dir",
        "db",
        "control_filename",
        "template_filename",
        "template_dir",
        "template_str",
        "control_settings",
        "save_runs_to_db",
        "save_runs_to_file",
        "lammpsin_key",
    ]

    def run_task(self, fw_spec):
        # Set directory for writing control file
        working_dir = fw_spec.get("working_dir", self.get("working_dir",
                                                          os.getcwd()))
        template_dir = None
        os.makedirs(working_dir, exist_ok=True)
        os.chdir(working_dir)

        db = self.get("db", None)
        save_to_db = self.get("save_runs_to_db", False)
        save_to_file = self.get("save_runs_to_file", True)

        # Set filename for control file
        control_filename = self.get("control_filename", "complex.lammpsin")

        # Set the settings parameter
        control_settings = self.get("control_settings", None)

        default_masses_list = fw_spec.get("default_masses", [])
        if default_masses_list:
            lammps_elements = lammps_mass_to_element(default_masses_list)
            if "X" not in lammps_elements:
                control_settings["dump_modify_elements"] = "element {}".format(
                    " ".join(lammps_elements)
                )

        # There are three different cases for input of the template:
        # Case 1: template as string
        if isinstance(self.get("template_str"), str):
            template_string = self.get("template_str")
            template_filename = None
        # Cases 2 and 3: template from file, with filename set by
        # template_type variable
        elif isinstance(self.get("template_filename"), str):
            template_filename = self.get("template_filename")
            # Case 2: template from file in /../templates/ directory
            if template_filename in TEMPLATE_TYPES:
                template_dir = TEMPLATE_DIR
            # Case 3: template from file provided by user
            elif isinstance(self.get("template_dir"), str):
                template_dir = self.get("template_dir")
            else:
                raise KeyError(
                    "Directory containing custom template file was not "
                    "specified; add as optional parameter"
                )
        else:
            raise KeyError(
                "Either template was not provided as a valid string or the "
                "path to a template text file was not provided."
            )

        if template_filename:
            template_path = os.path.join(template_dir, template_filename)
            with open(template_path) as file:
                template_string = file.read()

        write_lammps_inputs(
            working_dir,
            template_string,
            settings=control_settings,
            script_filename=control_filename,
        )

        smiles_list = fw_spec.get("smiles", [])
        n_mols_dict = fw_spec.get("nmols", {})
        num_mols_list = fw_spec.get("num_mols_list", [])
        lmp_box = fw_spec.get("box", None)
        num_atoms_per_mol_list = fw_spec.get("num_atoms_per_mol", [])
        recalc_masses_list = fw_spec.get("recalc_masses", [])

        run_doc = process_run(
            smiles_list,
            n_mols_dict,
            lmp_box,
            template_filename,
            control_settings
        )
        run_doc.update({"working_dir": working_dir})

        run_list = {}
        if save_to_db:
            run_db = get_db(input_db=db)
            run_id = run_db.insert_run(run_doc)
            run_list["lammps_run_id_list"] = run_id
            logger.info("Saved run info to db")

        if save_to_file:
            file = os.path.join(working_dir, "run_file.json")
            with open(file, "w") as run_file:
                run_file.write(json.dumps(run_doc, default=DATETIME_HANDLER))
            run_list["lammps_run_loc_list"] = file
            logger.info("Saved run info to json file")

        spec = {
            "smiles": smiles_list,
            "nmols": n_mols_dict,
            "num_mols_list": num_mols_list,
            "box": lmp_box,
            "num_atoms_per_mol": num_atoms_per_mol_list,
            "default_masses": default_masses_list,
            "recalc_masses": recalc_masses_list,
        }

        uid = self.get("lammpsin_key")
        set_dict = {f"lammpsin->{DEFAULT_KEY}": run_doc}
        if uid:
            set_dict[f"lammpsin->{uid}"] = run_doc
        mod_dict = {"_set": set_dict}
        if run_list:
            mod_dict.update({"_push": run_list})
        return FWAction(update_spec=spec, mod_spec=mod_dict, propagate=True)


# TODO: Write a firetask for writing a general data file
# TODO: Write a firetask for writing a general control file template \
#           (multiple fix options, general template)


@explicit_serialize
class WriteTleapScript(FiretaskBase):
    """
    Write a tleap script for generating a prmtop and inpcrd file for a molecule.
    Intended to be used to obtain the parameters from the GAFF force field.

    Args:
        working_dir (str): Path to the working directory. Defaults to current working
            directory.
        script_string (str): String containing the tleap script. If not provided, the
            ``template_filename`` and ``template_dir`` must be provided.
        template_filename (str): Name of the template file. Defaults to "gaff_tleap".
        template_dir (str): Path to the directory containing the template file. Defaults
            to the mispr/lammps/templates/ directory.
        tleap_settings (dict): Dictionary containing the settings for the tleap script.
            Used to update the ``TLEAP_SETTINGS`` dict in
            mispr/lammps/defaults.py file. Defaults to empty dict.
        script_filename (str): Name of the tleap script file. Defaults to "tleap.in".
    """
    _fw_name = "Write Tleap Script"
    required_params = []
    optional_params = [
        "working_dir",
        "script_string",
        "template_filename",
        "template_dir",
        "tleap_settings",
        "script_filename",
    ]

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

        tleap_settings = TLEAP_SETTINGS.update(
            self.get("tleap_settings", fw_spec.get("tleap_settings", {}))
        )

        script_string = Template(script_string).substitute(
            tleap_settings, **TLEAP_SETTINGS
        )

        script_filename = self.get("script_filename", "tleap.in")
        script_file_path = os.path.join(working_dir, script_filename)

        with open(script_file_path, "w") as script:
            script.write(script_string)


@explicit_serialize
class LabelFFDict(FiretaskBase):
    """
    Append molecule label to all atom labels in the force field dictionary. Ensures
    that no two molecules can share the same atom labels provided that each molecule
    has a unique label for the system.

    Args:
        mol (Molecule, optional): pymatgen Molecule object; if not provided, the
        ``prev_calc_molecule`` key from the ``fw_spec`` will be used.
        unlabeled_dict (dict, optional): Dictionary containing the force field parameters
            for the molecule; ignored if the ``ff_file`` is provided.
        ff_file (str, optional): Path to the json file containing the force field
            parameters for the molecule.
        working_dir (str, optional): Path to the working directory. Defaults to current
            working directory.
        label (str, optional): Label for the molecule. Defaults to the molecular formula.
        system_force_field_dict (dict, optional): Dictionary containing the force field
            parameters for the system that is used as the input for ``LammpsDataWrapper``
            object; intended to be obtained from the spec; if present in the spec,
            this will be ignored.
    """
    _fw_name = "Label FF Dict"
    required_params = []
    optional_params = [
        "mol",
        "unlabeled_dict",
        "ff_file",
        "working_dir",
        "label",
        "system_force_field_dict",
    ]

    def run_task(self, fw_spec):
        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        ff_file = self.get("ff_file")
        mol = self.get("mol", fw_spec.get("prev_calc_molecule"))

        if ff_file:
            with open(ff_file, "r") as file:
                unlabeled_dict = json.load(file)
                # mol = unlabeled_dict["Molecule"]
        else:
            if not isinstance(mol, Molecule):
                raise TypeError(
                    '"mol" input must be a PyMatGen Molecule object'
                )

            if not isinstance(self.get("unlabeled_dict"), dict):
                raise TypeError('"unlabeled_dict" input must be a dict')

            unlabeled_dict = self.get("unlabeled_dict")

        unlabeled_dict["Molecule"] = mol
        label = self.get("label", get_mol_formula(mol))
        labeled_dict = add_ff_labels_to_dict(unlabeled_dict, label)
        sys_ff_dict = fw_spec.get(
            "system_force_field_dict", self.get("system_force_field_dict", {})
        )
        sys_ff_dict[label] = labeled_dict
        return FWAction(update_spec={"system_force_field_dict": sys_ff_dict})


@explicit_serialize
class LabelFFDictFromDB(FiretaskBase):
    """
    Append molecule label to all atom labels in a force field dictionary obtained from
    the database. Ensures that no two molecules can share the same atom labels provided
    that each molecule has a unique label for the system.

    Args:
        filter (dict): Dictionary containing the filter used to query the database.
        working_dir (str, optional): Path to the working directory. Defaults to
            current working directory.
        db (dict or str, optional): Information for connecting to the database as
            either a dictionary of credentials or a path to a json file containing the
            credentials; if not provided, the database specified in the env will be used.
        label (str, optional): Label for the molecule. Defaults to the molecular formula.
        system_force_field_dict (dict, optional): Dictionary containing the force field
            parameters for the system that is used as the input for ``LammpsDataWrapper``
            object; intended to be obtained from the spec; if present in the spec, this
            will be ignored.
    """
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

        sys_ff_dict = fw_spec.get(
            "system_force_field_dict", self.get("system_force_field_dict", {})
        )

        sys_ff_dict[label] = labeled_ff_dict

        return FWAction(update_spec={"system_force_field_dict": sys_ff_dict})
