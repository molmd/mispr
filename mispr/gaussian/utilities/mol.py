# coding: utf-8


# Contains functions for processing molecules.

import os
import logging

from bson import ObjectId

from openbabel import pybel as pb
from openbabel import OBMolBondIter

from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.core.structure import Molecule

from mispr.gaussian.utilities.db_utilities import get_db

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)


def process_mol(operation_type, mol, local_opt=False, **kwargs):
    working_dir = kwargs.get("working_dir", os.getcwd())

    def get_db_():
        return get_db(kwargs["db"]) if "db" in kwargs else get_db()

    if operation_type == "get_from_mol":
        if not isinstance(mol, Molecule):
            raise Exception(
                "mol is not a Molecule object; either "
                "provide a Molecule object or use another "
                "operation type with its corresponding inputs"
            )
        output_mol = mol

    elif operation_type == "get_from_file" or operation_type == "get_from_gout_file":
        if not os.path.isabs(mol):
            file_path = os.path.join(working_dir, mol)
        else:
            file_path = mol
        if not os.path.exists(file_path):
            raise Exception(
                "mol is not a valid path; either provide a valid "
                "path or use another operation type with its "
                "corresponding inputs"
            )

        output_mol = Molecule.from_file(file_path)

    elif operation_type == "get_from_str":
        str_type = kwargs.get("str_type")
        if str_type is None:
            raise ValueError(
                "a mol string format must be specified to process " "the input string"
            )
        output_mol = Molecule.from_str(mol, str_type)

    elif operation_type == "get_from_mol_db":
        # mol = mol_inchi
        db = get_db_()
        mol_dict = db.retrieve_molecule(mol)
        if not mol_dict:
            raise Exception("mol is not found in the database")
        output_mol = Molecule.from_dict(mol_dict)

    elif operation_type == "get_from_gout":
        if not isinstance(mol, GaussianOutput):
            raise Exception(
                "mol is not a GaussianOutput object; either "
                "provide a GaussianOutput object or use another "
                "operation type with its corresponding inputs"
            )
        output_run = mol.as_dict()
        output_mol = Molecule.from_dict(output_run["output"]["molecule"])

    elif operation_type == "get_from_run_dict":
        if not isinstance(mol, dict) and "output" not in mol:
            raise Exception(
                "mol is not a GaussianOutput dictionary; either "
                "provide a GaussianOutput dictionary or use "
                "another operation type with its corresponding "
                "inputs"
            )
        output_mol = Molecule.from_dict(mol["output"]["output"]["molecule"])

    elif operation_type == "get_from_run_id":
        # mol = run_id
        db = get_db_()
        run = db.runs.find_one({"_id": ObjectId(mol)})
        if not run:
            raise Exception("Gaussian run is not in the database")
        mol_dict = run["output"]["output"]["molecule"]
        output_mol = Molecule.from_dict(mol_dict)

    elif operation_type == "get_from_run_query":
        # mol = {'inchi': inchi, 'type': type, 'functional': func,
        #        'basis': basis, 'phase': phase, ...}
        logger.info(
            "If the query criteria satisfy more than "
            "one document, the last updated one will "
            "be used. To perform a more specific "
            "search, provide the document id using "
            "gout_id"
        )
        db = get_db_()
        run = db.retrieve_run(**mol)
        if not run:
            raise Exception("Gaussian run is not in the database")
        run = max(run, key=lambda i: i["last_updated"])
        mol_dict = run["output"]["output"]["molecule"]
        output_mol = Molecule.from_dict(mol_dict)

    elif operation_type == "derive_molecule":
        # mol = {'operation_type': 'get_from_file', 'mol': file_path,
        #        'func_grp': func_group_name, ....}

        func_grp_name = mol.get("func_grp")
        if not func_grp_name:
            raise Exception(
                "No FG provided; Provide the name of the FG "
                "to be retrieved from the database"
            )
        db = get_db_()
        fg_dict = db.retrieve_fg(func_grp_name)
        if not fg_dict:
            raise Exception("FG is not found in the database")
        fg = Molecule(fg_dict["species"], fg_dict["coords"])

        output_mol = process_mol(
            operation_type=mol["operation_type"], mol=mol["mol"], **kwargs
        )
        output_mol.substitute(mol["index"], fg, mol["bond_order"])

    elif operation_type == "link_molecules":
        # TODO: add a checking step in the original molecule to make sure no
        #  overlapping happens
        # mol = {'operation_type': ['get_from_file', 'get_from_mol_db'],
        #        'mol': ['mol1.xyz', 'mol_inchi'],
        #        'index': [3, 5],
        #        'bond_order': 2}

        # mol = {'operation_type': ['get_from_file', 'derive_molecule'],
        #        'mol': ['mol2.xyz', {'operation_type':
        #        'get_from_mol_db, 'mol': inchi}],
        #        'index': [3, 5],
        #        'bond_order': 2}
        linking_mol = process_mol(
            operation_type=mol["operation_type"][0], mol=mol["mol"][0], **kwargs
        )
        linked_mol = process_mol(
            operation_type=mol["operation_type"][1], mol=mol["mol"][1], **kwargs
        )
        output_mol = linking_mol.link(
            linked_mol, mol["index"][0], mol["index"][1], mol["bond_order"]
        )

    else:
        raise ValueError(f"operation type {operation_type} is not supported")

    if local_opt:
        force_field = kwargs.get("force_field", "mmff94")
        steps = kwargs.get("steps", 500)
        output_mol = perform_local_opt(output_mol, force_field, steps)

    charge = kwargs.get("charge")
    charge = charge if charge is not None else output_mol.charge
    output_mol.set_charge_and_spin(charge)
    return output_mol


def perform_local_opt(mol, force_field="uff", steps=200):
    a = BabelMolAdaptor(mol)
    a.localopt(forcefield=force_field, steps=steps)
    mol_opt = a.pymatgen_mol
    return mol_opt


def label_atoms(mol):
    # supports all atom types and molecules as large as 999 atoms; does not
    # support clusters
    a = BabelMolAdaptor(mol)
    pm = pb.Molecule(a.openbabel_mol)
    mol_smiles = pm.write("smi").strip()
    mol_smiles_copy = mol_smiles.lower()
    counter = 0
    count_1 = ""
    count_2 = ""
    count_3 = ""
    smiles = mol_smiles.lower()
    atoms = [str(i).lower() for i in mol.species]
    sorted_atoms = sorted(atoms, key=lambda x: len(x), reverse=True)
    index_atom_map = {}
    h_index = [i for i, j in enumerate(mol_smiles_copy) if j == "h"]
    for atom in sorted_atoms:
        if atom == "h":
            smiles = smiles.replace(atom, "!", 1)
            continue
        index = smiles.find(atom)
        if index != -1:
            smiles = smiles.replace(atom, "!", 1)
            index_atom_map[index] = atom
    existing_atoms = [index_atom_map[i] for i in sorted(index_atom_map)]
    for ch in smiles:
        if ch == "!":
            flag = 0
            if len(count_1) in h_index:
                atom = "h"
                flag = 1
            else:
                atom = atoms[counter]
                while atom not in existing_atoms:
                    counter += 1
                    atom = atoms[counter]
            if atom.lower() == "h":
                i = 0
            else:
                i = counter
            if i < 10:
                count_1 += f"{i: >{len(atom)}}"
                count_2 += " " * len(atom)
                count_3 += " " * len(atom)
            elif i < 100:
                count_1 += f"{i // 10: >{len(atom)}}"
                count_2 += f"{i % 10: >{len(atom)}}"
                count_3 += " " * len(atom)
            else:
                count_1 += f"{i // 10 // 10: >{len(atom)}}"
                count_2 += f"{i % 100 // 10: >{len(atom)}}"
                count_3 += f"{i % 10: >{len(atom)}}"
            if flag == 0:
                counter += 1
        else:
            count_1 += " "
            count_2 += " "
            count_3 += " "
    print(f"{mol_smiles}\n{count_1}\n{count_2}\n{count_3}")


def get_bond_order_str(mol):
    """
    Finds bond order as a string by iterating over bonds of a molecule. First
    converts pymatgen mol to openbabel mol to use openbabel in finding bond
    order.
    """
    bond_order = {}
    a = BabelMolAdaptor(mol).openbabel_mol
    for bond in OBMolBondIter(a):
        atom_indices = tuple(sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]))
        order = bond.GetBondOrder()
        if order == 0:
            bond_order[atom_indices] = "U"
        elif order == 1:
            bond_order[atom_indices] = "S"
        elif order == 2:
            bond_order[atom_indices] = "D"
        elif order == 3:
            bond_order[atom_indices] = "T"
        elif order == 5:
            bond_order[atom_indices] = "A"
        else:
            raise TypeError("Bond order number {} is not understood".format(order))
    return bond_order
