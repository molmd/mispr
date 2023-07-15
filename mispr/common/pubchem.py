# coding: utf-8


# Implements a core class PubChemRunner for retrieving molecules from
# the PubChem database using a molecule name as a query criteria.

import os

import pubchempy as pcp

from pymatgen.core.structure import Molecule

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Aug 2022"
__version__ = "0.0.3"


class PubChemRunner:
    """
    Wrapper for retrieving molecules from PubChem database.
    """

    def __init__(self, abbreviation, working_dir=None):
        """
        Args:
            abbreviation (str): abbreviation to be used when saving molecule file.
            working_dir (str): working directory for saving the molecule file in;
                will use the current working directory if not specified
        """
        self.abbreviation = abbreviation
        self.working_dir = working_dir or os.getcwd()
        self.cid = None

    def get_mol(self, name, save_to_file=True, fmt="pdb", cleanup=True):
        """
        Wrapper function that searches for a molecule in the PubChem database,
        downloads it in the form of an SDF file, and converts the file
        to a pymatgen Molecule object.

        Args:
            name (str): name of the molecule to use for searching PubChem
            save_to_file (bool): whether to save the Molecule object to a file
            fmt (str): molecule file format if save_to_file is True
            cleanup (bool): whether to remove the intermediate sdf file

        Returns:
            Molecule: pymatgen Molecule object
        """
        self.download_sdf(name)
        molecule = self.convert_sdf_to_mol(save_to_file, fmt)
        if cleanup:
            self.cleanup()
        return molecule

    def download_sdf(self, name):
        """
        Downloads an SDF file from PubChem using a common name for the
        molecule as an identifier

        Args:
            name (str): name of the molecule to use for searching PubChem
        """
        cids = pcp.get_cids(name, record_type="3d")
        if cids:
            self.cid = cids[0]
            pcp.download(
                outformat="SDF",
                path=f"{self.working_dir}/{self.abbreviation}_{self.cid}.sdf",
                identifier=self.cid,
                record_type="3d",
                overwrite=True,
            )
        else:
            raise ValueError("No matching molecule found in the PubChem database")

    def convert_sdf_to_mol(self, save_to_file, fmt):
        """
        Converts an SDF file to a pymatgen Molecule object

        Args:
            save_to_file (bool): whether to save the Molecule object to a file
            fmt (str): molecule file format if save_to_file is True

        Returns:
            Molecule: pymatgen Molecule object
        """
        mol = Molecule.from_file(
            f"{self.working_dir}/{self.abbreviation}_{self.cid}.sdf"
        )
        if save_to_file:
            mol.to(fmt, f"{self.working_dir}/{self.abbreviation}_{self.cid}.{fmt}")
        return mol

    def cleanup(self):
        """
        Deletes the sdf file downloaded from PubChem
        """
        os.remove(f"{self.working_dir}/{self.abbreviation}_{self.cid}.sdf")
