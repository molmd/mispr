# coding: utf-8


# Contains functions for creating db schema.

import logging

from openbabel import pybel as pb

from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.core.structure import Molecule

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)


def get_chem_schema(mol):
    mol_dict = mol.as_dict()
    comp = mol.composition
    a = BabelMolAdaptor(mol)
    pm = pb.Molecule(a.openbabel_mol)
    # svg = pm.write('svg')
    mol_dict.update(
        {
            "smiles": pm.write("smi").strip(),
            "inchi": pm.write("inchi").strip("\n"),
            "formula": comp.formula,
            "formula_pretty": comp.reduced_formula,
            "formula_anonymous": comp.anonymized_formula,
            "formula_alphabetical": comp.alphabetical_formula,
            "chemsys": comp.chemical_system,
            "nsites": mol.num_sites,
            "nelements": len(comp.chemical_system.replace("-", " ").split(" ")),
            "is_ordered": mol.is_ordered,
            "is_valid": mol.is_valid(),
        }
    )
    return mol_dict


def get_mol_formula(mol):
    mol_schema = get_chem_schema(mol)
    return mol_schema["formula_alphabetical"].replace(" ", "")


def get_job_name(mol, name):
    if not isinstance(mol, Molecule):
        job_name = "{}_{}".format(mol, name)
    else:
        mol_formula = get_mol_formula(mol)
        job_name = "{}_{}".format(mol_formula, name)
    return job_name
