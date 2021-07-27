# coding: utf-8


# Contains functions for processing rdkit molecules.

import os
import random
import logging

from mispr.gaussian.utilities.mol import get_bond_order_str

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)


def get_rdkit_mol(mol, sanitize=True, remove_h=False):
    """
    Converts a pymatgen mol object to RDKit rdmol object. Uses RDKit to perform
    the conversion <http://rdkit.org>. Accounts for aromaticity.
    """
    try:
        import rdkit
        from rdkit import Chem
        from rdkit.Geometry import Point3D
    except ModuleNotFoundError:
        raise ImportError("This function requires RDKit to be installed.")

    mol_species = [str(s) for s in mol.species]
    mol_coords = mol.cart_coords

    rdkit_mol = Chem.rdchem.EditableMol(Chem.rdchem.Mol())
    conformer = Chem.Conformer(len(mol_species))

    for index, (specie, coord) in enumerate(zip(mol_species, mol_coords)):
        rdkit_atom = Chem.rdchem.Atom(specie)
        rdkit_mol.AddAtom(rdkit_atom)
        conformer.SetAtomPosition(index, Point3D(*coord))

    rdkit_bonds = Chem.rdchem.BondType
    bond_order_mapping = {
        "U": rdkit_bonds.UNSPECIFIED,
        "S": rdkit_bonds.SINGLE,
        "D": rdkit_bonds.DOUBLE,
        "T": rdkit_bonds.TRIPLE,
        "A": rdkit_bonds.AROMATIC,
    }
    bond_orders = get_bond_order_str(mol)

    for bond, bond_order in bond_orders.items():
        order = bond_order_mapping[bond_order]
        rdkit_mol.AddBond(bond[0] - 1, bond[1] - 1, order)

    rdkit_mol = rdkit_mol.GetMol()
    if sanitize:
        # if sanitization fails, inform the user and proceed normally
        try:
            Chem.SanitizeMol(rdkit_mol)
        except Exception as e:
            logger.info("Could not sanitize rdkit mol: {}".format(e))
    rdkit_mol.AddConformer(conformer, assignId=False)
    if remove_h:
        rdkit_mol = Chem.RemoveHs(rdkit_mol, sanitize=sanitize)

    return rdkit_mol


def calc_energy(rdkit_mol, maxIters=200):
    from rdkit.Chem import AllChem

    AllChem.UFFOptimizeMolecule(rdkit_mol, maxIters)
    ff = AllChem.UFFGetMoleculeForceField(rdkit_mol)
    ff.Minimize(maxIts=maxIters)
    E = ff.CalcEnergy()
    return E


def draw_rdkit_mol(rdkit_mol, filename="mol.png", working_dir=None):
    from rdkit.Chem import Draw, AllChem

    working_dir = working_dir or os.getcwd()
    AllChem.Compute2DCoords(rdkit_mol)
    Draw.MolToFile(rdkit_mol, os.path.join(working_dir, filename))


def draw_rdkit_mol_with_highlighted_bonds(
    rdkit_mol, bonds, filename="mol.png", colors=None, working_dir=None
):
    def _generate_color():
        bond_color = (random.random(), random.random(), random.random())
        # prevent the generation of a white color
        if bond_color == (1.0, 1.0, 1.0):
            bond_color = _generate_color()
        return bond_color

    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    working_dir = working_dir or os.getcwd()
    rdDepictor.SetPreferCoordGen(True)
    d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    rdkit_mol_copy = Chem.Mol(rdkit_mol.ToBinary())
    # optimize 2D depiction
    rdDepictor.Compute2DCoords(rdkit_mol_copy, bondLength=1.5)
    d.SetFontSize(0.6 * d.FontSize())  # reduce font size of bond number

    highlighted_bonds = []
    for i, bond in enumerate(bonds):
        ind = rdkit_mol_copy.GetBondBetweenAtoms(*bond).GetIdx()
        rdkit_mol_copy.GetBondWithIdx(ind).SetProp("bondNote", str(i))
        highlighted_bonds.append(ind)

    # use randomly generated colors if no or wrong number of colors are provided
    if not colors or len(colors) < len(bonds):
        colors = []
        for i in range(len(bonds)):
            colors.append(_generate_color())

    bond_colors = {}
    for i, bond in enumerate(highlighted_bonds):
        bond_colors[bond] = colors[i]

    rdMolDraw2D.PrepareAndDrawMolecule(
        d,
        rdkit_mol_copy,
        highlightBonds=highlighted_bonds,
        highlightBondColors=bond_colors,
    )

    d.DrawMolecule(rdkit_mol_copy)
    d.FinishDrawing()
    d.WriteDrawingText(os.path.join(working_dir, filename))
    return colors
