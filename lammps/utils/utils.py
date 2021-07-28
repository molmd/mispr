import os
import numpy as np
from collections import OrderedDict
import infrastructure.lammps.firetasks.run as ilfr
import infrastructure.gaussian.utils.utils as iguu

def add_ff_labels_to_BADI_lists(list, label):
    """
    Adds extra string to the end of all atom type labels in lists containing
    information about Bonds, Angles, Dihedrals, or Impropers (BADI). This
    function is intended to be used through the add_ff_labels_to_dict()
    function.
    :param list: [List] The value from ff_dict using one of the following keys:
                 'Bonds', 'Angles', 'Dihedrals', or 'Impropers'.
                 The form of this list should be as follows:
                 [{'coeffs': [Float, ...], 'types': [(Str, ...), ...]}, ...]
    :param label: [Str] A label for the molecular species that is unique for
                  the system being created.
    :return:
    """
    output_badi_list = []
    for dict in list:
        new_types = []
        for type in dict['types']:
            new_types.append(tuple(atom + label for atom in type))
        output_badi_list.append({'coeffs': dict['coeffs'], 'types': new_types})
    return output_badi_list


def add_ff_labels_to_dict(ff_dict, label):
    """

    :param ff_dict:
                    {'Molecule':   pmg.Molecule,
                     'Labels':     List,
                     'Masses':     OrderedDict,
                     'Nonbond':    List,
                     'Bonds':      [{'coeffs': [a, b], 'types': [('x1', 'x2'), ...]}, ...],
                     'Angles':     [{'coeffs': [a, b], 'types': [('x1', 'x2', 'x3'), ...]}, ...],
                     'Dihedrals':  [{'coeffs': [a, b, c], 'types': [('x1', 'x2', 'x3', 'x4), ...]}, ...],
                     'Impropers':  [{'coeffs': [a, b, c], 'types': [('x1', 'x2', 'x3', 'x4), ...]}, ...],
                     'Improper Topologies': List,
                     'Charges':    np.Array,
                     ...}
    :param label:
    :return:
    """

    output_labels = [old_label + label for old_label in ff_dict['Labels']]

    output_masses = OrderedDict()
    for atom_type, mass in ff_dict['Masses'].items():
        output_masses[atom_type + label] = mass

    output_bonds = add_ff_labels_to_BADI_lists(ff_dict['Bonds'], label)
    output_angles = add_ff_labels_to_BADI_lists(ff_dict['Angles'], label)
    output_dihedrals = add_ff_labels_to_BADI_lists(ff_dict['Dihedrals'], label)
    output_impropers = add_ff_labels_to_BADI_lists(ff_dict['Impropers'], label)

    output_ff_dict = {'Molecule':   ff_dict['Molecule'],
                      'Labels':     output_labels,
                      'Masses':     output_masses,
                      'Nonbond':    ff_dict['Nonbond'],
                      'Bonds':      output_bonds,
                      'Angles':     output_angles,
                      'Dihedrals':  output_dihedrals,
                      'Impropers':  output_impropers,
                      'Improper Topologies': ff_dict['Improper Topologies'],
                      'Charges':    ff_dict['Charges']}

    return output_ff_dict


def get_db(input_db=None):
    from infrastructure.lammps.database import LammpsSysDb
    if not input_db:
        input_db = f"{ilfr.CONFIG_PATH}/db.json"
        if not os.path.isfile(input_db):
            raise FileNotFoundError(
                "Please provide the database configurations")
    if isinstance(input_db, dict):
        db = LammpsSysDb(**input_db)
    # else:
    #     db = LammpsSysDb.from_db_file(input_db)

    return db


def process_ff_doc(parameter_dict, method=None, doi=None, **kwargs):
    mol = parameter_dict.pop("Molecule")
    ff_dict = iguu.get_chem_schema(mol)
    ff_dict.update(parameter_dict)
    ff_dict["method"] = method
    ff_dict["doi"] = doi
    ff_dict.update(kwargs)
    return ff_dict


def process_run(smiles, nmols, box, template_filename, control_settings):
    if box is not None:
        box_setting = box.as_dict()
    else:
        box_setting = {}
    run_dict = {"smiles": smiles,
                "nmols": nmols,
                "box": box_setting,
                "job_type": template_filename}
    return run_dict

if __name__ == "__main__":
    Spce_label = ''
    test_mass = OrderedDict({'ow' + Spce_label: 16.000, 'hw' + Spce_label: 1.008})
    test_mass_2 = OrderedDict()

    for key in test_mass.keys():
        test_mass_2[key + 'H2O'] = test_mass[key]

    print(test_mass)
    print(test_mass_2, '\n')

    test_bonds = [{'coeffs': [553.0, 1], 'types': [('ow' + Spce_label, 'hw' + Spce_label)]},
                  {'coeffs': [552.0, 1], 'types': [('ow1' + Spce_label, 'hw1' + Spce_label)]}]

    test_bonds_2 = add_ff_labels_to_BADI_lists(test_bonds, 'H2O')

    print(test_bonds)
    print(test_bonds_2, '\n')

    Spce_param_dict = {'Molecule': 'Spce_molecule',
                       'Labels': ['ow' + Spce_label, 'hw' + Spce_label, 'hw' + Spce_label],
                       'Masses': OrderedDict({'ow' + Spce_label: 16.000, 'hw' + Spce_label: 1.008}),
                       'Nonbond': [[0.155394259, 3.16555789], [0.0, 0.0]],
                       'Bonds': [{'coeffs': [553.0, 1], 'types': [('ow' + Spce_label, 'hw' + Spce_label)]}],
                       'Angles': [{'coeffs': [100.0, 109.47],
                                   'types': [('hw' + Spce_label, 'ow' + Spce_label, 'hw' + Spce_label)]}],
                       'Dihedrals': [],
                       'Impropers': [],
                       'Improper Topologies': None,
                       'Charges': np.asarray([-0.8476, 0.4238, 0.4238])}
    Spce_param_dict_2 = add_ff_labels_to_dict(Spce_param_dict, 'H2O')

    print(Spce_param_dict)
    print(Spce_param_dict_2)