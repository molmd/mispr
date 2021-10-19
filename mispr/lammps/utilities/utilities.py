# coding: utf-8


# Defines lammps utility functions.

import os

from collections import OrderedDict

from fireworks.fw_config import CONFIG_FILE_DIR

from mispr.lammps.firetasks.run import RunLammpsFake
from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.gaussian.utilities.fw_utilities import get_list_fireworks_and_tasks

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"


def add_ff_labels_to_BADI_lists(ff_list, label):
    """
    Adds extra string to the end of all atom type labels in lists containing
    information about Bonds, Angles, Dihedrals, or Impropers (BADI). This
    function is intended to be used through the add_ff_labels_to_dict()
    function.
    :param ff_list: [List] The value from ff_dict using one of the following keys:
                 'Bonds', 'Angles', 'Dihedrals', or 'Impropers'.
                 The form of this list should be as follows:
                 [{'coeffs': [Float, ...], 'types': [(Str, ...), ...]}, ...]
    :param label: [Str] A label for the molecular species that is unique for
                  the system being created.
    :return:
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
    :param ff_dict:
        {'Molecule':
            pmg.Molecule,
            'Labels': List,
            'Masses': OrderedDict,
            'Nonbond': List,
            'Bonds': [{'coeffs': [a, b], 'types': [('x1', 'x2'), ...]}, ...],
            'Angles': [{'coeffs': [a, b], 'types': [('x1', 'x2', 'x3'), ...]}, ...],
            'Dihedrals': [{'coeffs': [a, b, c], 'types': [('x1', 'x2', 'x3', 'x4), ...]}, ...],
            'Impropers': [{'coeffs': [a, b, c], 'types': [('x1', 'x2', 'x3', 'x4), ...]}, ...],
            'Improper Topologies': List,
            'Charges': np.Array,
            ...}
    :param label:
    :return:
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
