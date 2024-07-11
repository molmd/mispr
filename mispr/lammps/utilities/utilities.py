"""Define lammps utility functions."""

import os
import json
import math
import numpy as np
import pandas as pd

from collections import OrderedDict

from fireworks.fw_config import CONFIG_FILE_DIR

from pymatgen.core.structure import Molecule
from pymatgen.io.lammps.data import Topology

from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.gaussian.utilities.fw_utilities import get_list_fireworks_and_tasks

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.4"


def add_ff_labels_to_BADI_lists(ff_list, label):
    """
    Add extra string to the end of all atom type labels in lists containing information
    about Bonds, Angles, Dihedrals, or Impropers (BADI). This function is intended to be
    used through the ``add_ff_labels_to_dict``.

    Args:
        ff_list (List): The value from ``ff_dict`` using one of the following keys:
            'Bonds', 'Angles', 'Dihedrals', or 'Impropers'. The form of this list
            should be as follows:

            .. code-block:: python

                [{'coeffs': [Float, ...], 'types': [(Str, ...), ...]}, ...]

        label (str): A label for the molecular species that is unique for the system
            being created.
    """
    output_badi_list = []
    for dict_ in ff_list:
        new_types = []
        for type in dict_["types"]:
            new_types.append(tuple(atom + label for atom in type))
        output_badi_list.append({"coeffs": dict_["coeffs"], "types": new_types})
    return output_badi_list


def add_ff_labels_to_dict(ff_dict, label):
    """
    Args:
        ff_dict: A dictionary containing the force field information for a molecule.
            The dictionary should have the following form:

            .. code-block:: python

                {
                    "Molecule": pmg.Molecule,
                    "Labels": List,
                    "Masses": OrderedDict,
                    "Nonbond": List,
                    "Bonds": [{'coeffs': [a, b], 'types': [('x1', 'x2'), ...]}, ...],
                    "Angles": [{'coeffs': [a, b], 'types': [('x1', 'x2', 'x3'), ...]}, ...],
                    "Dihedrals": [{'coeffs': [a, b, c], 'types': [('x1', 'x2', 'x3', 'x4), ...]}, ...],
                    "Impropers": [{'coeffs': [a, b, c], 'types': [('x1', 'x2', 'x3', 'x4), ...]}, ...],
                    "Improper Topologies": List,
                    "Charges": np.Array,
                    ...
                }

        label (str):
    """
    output_labels = [old_label + label for old_label in ff_dict["Labels"]]

    output_masses = OrderedDict()
    for atom_type, mass in ff_dict["Masses"].items():
        output_masses[atom_type + label] = mass

    output_bonds = add_ff_labels_to_BADI_lists(ff_dict["Bonds"], label)
    output_angles = add_ff_labels_to_BADI_lists(ff_dict["Angles"], label)
    output_dihedrals = add_ff_labels_to_BADI_lists(ff_dict["Dihedrals"], label)
    output_impropers = add_ff_labels_to_BADI_lists(ff_dict["Impropers"], label)

    output_ff_dict = {
        "Molecule": ff_dict["Molecule"],
        "Labels": output_labels,
        "Masses": output_masses,
        "Nonbond": ff_dict["Nonbond"],
        "Bonds": output_bonds,
        "Angles": output_angles,
        "Dihedrals": output_dihedrals,
        "Impropers": output_impropers,
        "Improper Topologies": ff_dict["Improper Topologies"],
        "Charges": ff_dict["Charges"],
    }
    return output_ff_dict


def get_db(input_db=None):
    from mispr.lammps.database import LammpsSysDb

    if not input_db:
        input_db = f"{CONFIG_FILE_DIR}/db.json"
        if not os.path.isfile(input_db):
            raise FileNotFoundError("Please provide the database configurations")
    if isinstance(input_db, dict):
        db = LammpsSysDb(**input_db)
    else:
        db = LammpsSysDb.from_db_file(input_db)
    return db


def process_ff_doc(parameter_dict, method=None, doi=None, **kwargs):
    mol = parameter_dict.pop("Molecule")
    ff_dict = get_chem_schema(mol)
    ff_dict.update(parameter_dict)
    ff_dict["method"] = method
    ff_dict["doi"] = doi
    ff_dict.update(kwargs)
    return ff_dict


def process_run(smiles, nmols, box, template_filename, control_settings=None):
    if box is not None:
        box_setting = box.as_dict()
    else:
        box_setting = {}
    run_dict = {
        "smiles": smiles,
        "nmols": nmols,
        "box": box_setting,
        "job_type": template_filename,
        "control_settings": control_settings,
    }
    return run_dict


def run_fake_lammps(workflow, ref_dirs, control_filenames=None):
    from mispr.lammps.firetasks.run import RunLammpsFake

    list_fireworks_and_tasks = get_list_fireworks_and_tasks(
        workflow, task_substring=["RunLammps"]
    )
    if not control_filenames:
        control_filenames = ["complex.lammpsin"] * len(ref_dirs)

    for ind, (i_firework, i_task) in enumerate(list_fireworks_and_tasks):
        workflow.fws[i_firework].tasks[i_task] = RunLammpsFake(
            ref_dir=ref_dirs[ind], control_filename=control_filenames[ind]
        )
    return workflow


def lammps_mass_to_element(lammps_masses):
    """
    Create a dict for mapping atom mass to element.

    Args:
        lammps_masses (list): List of masses in lammps units.

    Returns:
        dict
    """
    with open(os.path.join(os.path.dirname(__file__), "../data/masses.json")) as f:
        masses = json.load(f)

    elements = ["X"] * len(lammps_masses)
    for ind, mass in enumerate(lammps_masses):
        for item in masses.items():
            if math.isclose(mass, item[1], abs_tol=0.01):
                elements[ind] = item[0]
    return elements

def single_lammps_mass_to_element(lammps_mass):
    """
    Convert single lammps mass to element.

    Args:
        lammps_mass (float): Mass in lammps units.

    Returns:
        str
    """
    with open(os.path.join(os.path.dirname(__file__), "../data/masses.json")) as f:
        masses = json.load(f)

    for item in masses.items():
        if math.isclose(lammps_mass, item[1], abs_tol=0.01):
            element = item[0]
    return element


def unwrap_dump(dump_object, mol_file_path_dict, num_mol_dict, sorted_mol_names):
    """
    Unwraps molecules in a dump object. Sometimes, molecules can be
    split across periodic boundaries in a dump file. This function
    will unwrap the molecules so that they are not split across
    periodic boundaries.

    Args:
        dump_object (Dump): A pymatgen ``Dump`` object that will be
            unwrapped.
        mol_file_path_dict (dict): A dictionary with the molecule name
            as the key and the file path to the molecule in string form
            as the value.
        num_mol_dict (dict): A dictionary with the molecule name as the
            key and the number of molecules in the dump file in integer
            form as the value.
        sorted_mol_names (list): A list of the molecule names sorted in
            the order that they appear in the dump file. This should
            match the order of the molecules as they appear in the data
            file used to run the LAMMPS simulation.
    """
    dump_object.data = dump_object.data.sort_values(by="id")
    octant_bounds = np.mean(dump_object.box.bounds, axis=1)
    octant_lengths = np.diff(dump_object.box.bounds, axis=1)
    dump_object.data["oct_x"] = dump_object.data["x"].apply(lambda x: "h" if x > octant_bounds[0] else "l")
    dump_object.data["oct_y"] = dump_object.data["y"].apply(lambda y: "h" if y > octant_bounds[1] else "l")
    dump_object.data["oct_z"] = dump_object.data["z"].apply(lambda z: "h" if z > octant_bounds[2] else "l")
    dump_object.data["octant"] = dump_object.data["oct_x"] + dump_object.data["oct_y"] + dump_object.data["oct_z"]

    cum_num_mol_list = list(np.cumsum([0] + [num_mol_dict[mol_name] for mol_name in sorted_mol_names]) + 1)
    print(cum_num_mol_list)

    for i, mol_name in enumerate(sorted_mol_names):
        cur_mol = Molecule.from_file(mol_file_path_dict[mol_name])
        cur_bonds = Topology.from_bonding(cur_mol).topologies["Bonds"]
        for j in range(cum_num_mol_list[i], cum_num_mol_list[i + 1]):
            cur_mol_df = dump_object.data[dump_object.data["mol"] == j]
            cur_main_octant = cur_mol_df["octant"].mode().values[0]
            for rep in range(len(cur_bonds)):
                if np.all(cur_mol_df["octant"] == cur_main_octant):
                    break
                else:
                    num_broken_bonds = 0
                    for k, bond in enumerate(cur_bonds):
                        cur_bond_atoms = cur_mol_df.iloc[bond]
                        print(cur_bond_atoms)
                        print(bond)
                        cur_bond_length = np.linalg.norm(cur_bond_atoms.iloc[0][["x", "y", "z"]] - cur_bond_atoms.iloc[1][["x", "y", "z"]])
                        print(cur_bond_length)
                        if abs(cur_bond_length) > abs(octant_bounds.max()):
                            num_broken_bonds += 1
                        for l, coord in enumerate(["x", "y", "z"]):
                            cur_coord_diff = cur_bond_atoms[coord].diff().dropna().values[0]
                            print(cur_coord_diff)
                            if abs(cur_coord_diff) > octant_bounds[l]:
                                # print(cur_bond_atoms)
                                # print(bond)
                                atom_ids_to_change = cur_bond_atoms[cur_bond_atoms[f"oct_{coord}"] != cur_main_octant[l]].index
                                print(cur_main_octant[l])
                                print(atom_ids_to_change)
                                if cur_main_octant[l] == "h":
                                    print("h")
                                    print(dump_object.data.loc[atom_ids_to_change])
                                    dump_object.data.loc[atom_ids_to_change, coord] = cur_bond_atoms.loc[atom_ids_to_change, coord] + octant_lengths[l]
                                    cur_mol_df.loc[atom_ids_to_change, coord] = cur_bond_atoms.loc[atom_ids_to_change, coord] + octant_lengths[l]
                                    cur_mol_df.loc[atom_ids_to_change, f"oct_{coord}"] = "h"
                                    cur_mol_df.loc[atom_ids_to_change, "octant"] = cur_bond_atoms["oct_x"] + cur_bond_atoms["oct_y"] + cur_bond_atoms["oct_z"]
                                    print(dump_object.data.loc[atom_ids_to_change])
                                    print(cur_mol_df.loc[atom_ids_to_change])
                                elif cur_main_octant[l] == "l":
                                    print("l")
                                    print(dump_object.data.loc[atom_ids_to_change])
                                    dump_object.data.loc[atom_ids_to_change, coord] = cur_bond_atoms.loc[atom_ids_to_change, coord] - octant_lengths[l]
                                    cur_mol_df.loc[atom_ids_to_change, coord] = cur_bond_atoms.loc[atom_ids_to_change, coord] - octant_lengths[l]
                                    cur_mol_df.loc[atom_ids_to_change, f"oct_{coord}"] = "l"
                                    cur_mol_df.loc[atom_ids_to_change, "octant"] = cur_bond_atoms["oct_x"] + cur_bond_atoms["oct_y"] + cur_bond_atoms["oct_z"]
                                    print(dump_object.data.loc[atom_ids_to_change])
                                    print(cur_mol_df.loc[atom_ids_to_change])
                    if num_broken_bonds == 0:
                        break


def concat_slab(slab_file,
                bulk_file,
                output_file="interface.xyz",
                slab_location="bottom",
                origin=[0.0, 0.0, 0.0],
                p_x=0.0,
                p_y=0.0,
                interfacial_separation=2.0,
                output_first_atoms="bulk",
                flip_slab=False,
                working_dir=None):
    """
    Concatenate the slab and bulk xyz files to form an interface xyz 
    file. The slab is assumed to be normal to the z direction and can be
    added to the top or the bottom of the bulk. This function can be run
    more than once to create a sandwich structure of slab-bulk-slab.

    Args:
        slab_file (str): The path to the slab xyz file.
        bulk_file (str): The path to the bulk xyz file.
        output_file (str, optional): The path to the output xyz file.
            Defaults to "interface.xyz".
        slab_location (str, optional): The location of the slab relative
            to the z axis. Either "bottom" or "top" based on whether the
            slab will have the lowest z coords or the highest z coords.
            Defaults to "bottom".
        origin (list, optional): The origin of the interfacial system in
            the output file. Defaults to [0.0, 0.0, 0.0].
        p_x (float, optional): The extra space needed to establish the
            periodic boundary in the x direction. For simple fcc
            crystals, this is the atomic radius. Defaults to 0.0.
        p_y (float, optional): The same as p_x except in the y
            direction.
        interfacial_separation (float, optional): The separation between
            the slab and the bulk. Defaults to 2.0.
        output_first_atoms (str): The first atoms to be written to the
            output file. Either "slab" or "bulk".
        flip_slab (bool, optional): Whether to flip the slab in the z
            direction.
        working_dir (str, optional): The working directory where the
            output file will be written. Defaults to None, in which case
            the current working directory is used.
        
    Returns:
        str: The path to the output xyz file.

    Raises:
        ValueError: If slab_location is not "bottom" or "top".
        ValueError: If output_first_atoms is not "slab" or "bulk".
    """

    # Read the slab and bulk xyz files
    slab_df = pd.read_csv(
        slab_file, 
        delim_whitespace=True, 
        skiprows=2, 
        names=["element", "x", "y", "z"]
    )
    bulk_df = pd.read_csv(
        bulk_file, 
        delim_whitespace=True, 
        skiprows=2, 
        names=["element", "x", "y", "z"]
    )

    if slab_location not in ["bottom", "top"]:
        raise ValueError("slab_location must be either 'bottom' or 'top'")
    
    if output_first_atoms not in ["slab", "bulk"]:
            raise ValueError("output_first_atoms must be either 'slab' or 'bulk'")

    # Reorient the slab to start at the origin; will only reorient the z
    # coordinates of the slab to the origin if it is at the bottom.
    if flip_slab:
        slab_df["z"] *= -1

    slab_x_min = slab_df["x"].min()
    slab_y_min = slab_df["y"].min()
    slab_z_min = slab_df["z"].min()
    slab_x_max = slab_df["x"].max()
    slab_y_max = slab_df["y"].max()
    slab_z_max = slab_df["z"].max()

    slab_df["x"] += origin[0] - slab_x_min
    slab_df["y"] += origin[1] - slab_y_min
    if slab_location == "bottom":
        slab_df["z"] += origin[2] - slab_z_min

    # Reorient the bulk to start at the top of the slab, including the 
    # interfacial separation; only does so for the z coordinates if the
    # slab is on bottom.
    bulk_x_min = bulk_df["x"].min()
    bulk_y_min = bulk_df["y"].min()
    bulk_z_min = bulk_df["z"].min()
    bulk_x_max = bulk_df["x"].max()
    bulk_y_max = bulk_df["y"].max()
    bulk_z_max = bulk_df["z"].max()

    bulk_df["x"] += origin[0] - bulk_x_min
    bulk_df["y"] += origin[1] - bulk_y_min
    if slab_location == "bottom":
        bulk_df["z"] += origin[2] + slab_z_max - slab_z_min + interfacial_separation - bulk_z_min

    # Rescale the x and y coordinates of the bulk to account for the
    # periodic boundary conditions.
    bulk_df["x"] *= (slab_x_max - slab_x_min + p_x) / (bulk_x_max - bulk_x_min)
    bulk_df["y"] *= (slab_y_max - slab_y_min + p_y) / (bulk_y_max - bulk_y_min)

    # If the slab is on top, then assign the z coordinates of the slab
    # and bulk.
    if slab_location == "top":
        bulk_df["z"] += origin[2] - bulk_z_min
        slab_df["z"] += origin[2] + bulk_z_max - bulk_z_min + interfacial_separation - slab_z_min

    # Concatenate the slab and bulk dataframes.
    if output_first_atoms == "slab":
        interface_df = pd.concat([slab_df, bulk_df])
    elif output_first_atoms == "bulk":
        interface_df = pd.concat([bulk_df, slab_df])
    else:
        raise ValueError("output_first_atoms must be either 'slab' or 'bulk'")
    
    print(interface_df)
    
    # Write the interface dataframe to the output xyz file.
    if working_dir is None:
        working_dir = os.getcwd()
    output_path = os.path.join(working_dir, output_file)
    natoms = interface_df.shape[0]
    with open(output_path, "w") as f:
        f.write(f"{natoms}\n")
        f.write(f"interface from slab {slab_file} and bulk {bulk_file}\n")
        interface_df.to_csv(f, sep=" ", header=False, index=False)

    return output_path