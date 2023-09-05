# coding: utf-8


# Defines firetasks for writing Gaussian input files.

import os

from copy import deepcopy

from pymatgen.io.gaussian import GaussianInput
from pymatgen.core.structure import IMolecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.4"


@explicit_serialize
class WriteInput(FiretaskBase):
    """
    Write Gaussian input file for a molecule/cluster in the current
    working directory.

    Optional Args:
        gaussian_input (GaussianInput): a GaussianInput object to be
            written to a file; can either be provided directly or
            passed through fw_spec; if both are provided, will take
            the one provided directly
        molecule (Molecule): a Molecule object to be used in creating
            Gaussian input file; can be provided as an optional parameter
            or passed through fw_spec; if both are provided, will check
            if their graphs are isomorphic; if they are isomorphic,
            will use the molecule passed through fw_spec, otherwise
            will take the user input; note that molecule is ignored if
            a gaussian_input is available
        gaussian_input_params (dict): a dictionary of parameters to be
            used in creating the Gaussian input file; ignored if a
            gaussian_input is provided
        input_file (str): name of the input file to be written; default
            is mol.com
        cart_coords (bool): whether to write cartesian coordinates or
            not; default is True
        oxidation_states (dict): a dictionary of element symbols and
            their oxidation states used in setting the charge on the
            molecule
    """

    _fw_name = "Write Gaussian Input File"

    required_params = []
    optional_params = [
        "gaussian_input",
        "molecule",
        "gaussian_input_params",
        "input_file",
        "cart_coords",
        "oxidation_states",
    ]

    def _update_charge(self, mol):
        """
        Calculates the charge of a molecule/cluster using the oxidation state
        of individual elements.
        """
        if self.get("oxidation_states"):
            mol_copy = deepcopy(mol)
            mol_copy.add_oxidation_state_by_element(self.get("oxidation_states"))
            mol_copy.set_charge_and_spin(super(IMolecule, mol_copy).charge)
            self["gaussian_input_params"] = {
                **self.get("gaussian_input_params", {}),
                "charge": int(mol_copy.charge),
            }

    def run_task(self, fw_spec):
        working_dir = os.getcwd()

        input_file = self.get("input_file", "mol.com")
        input_path = os.path.join(working_dir, input_file)

        # if a full Gaussian object is provided
        if self.get("gaussian_input") and isinstance(
            self.get("gaussian_input"), GaussianInput
        ):
            gaussin = self["gaussian_input"]

        # if a full Gaussian object is being passed through fw_spec
        elif fw_spec.get("gaussian_input") and isinstance(
            fw_spec.get("gaussian_input"), GaussianInput
        ):
            gaussin = fw_spec.get("gaussian_input")

        # if a molecule is being passed through fw_spec
        elif fw_spec.get("prev_calc_molecule"):
            prev_calc_mol = fw_spec.get("prev_calc_molecule")
            # if a molecule is also passed as an optional parameter
            if self.get("molecule"):
                mol = self.get("molecule")
                mol_graph = MoleculeGraph.with_local_env_strategy(
                    mol, OpenBabelNN(), reorder=False, extend_structure=False
                )
                prev_mol_graph = MoleculeGraph.with_local_env_strategy(
                    prev_calc_mol, OpenBabelNN(), reorder=False, extend_structure=False
                )
                if mol_graph.isomorphic_to(prev_mol_graph):
                    mol = prev_calc_mol
                else:
                    print(
                        "Not using prev_calc_mol as it is not isomorphic to passed molecule!"
                    )
            else:
                mol = prev_calc_mol
            self._update_charge(mol)
            gaussin = GaussianInput(mol, **self.get("gaussian_input_params", {}))
        # if a molecule is only included as an optional parameter
        elif self.get("molecule"):
            self._update_charge(self.get("molecule"))
            gaussin = GaussianInput(
                self.get("molecule"), **self.get("gaussian_input_params", {})
            )
        # if no molecule is present raise an error
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )
        gaussin.write_file(input_path, self.get("cart_coords", True))
        # delete fw_spec to avoid pymatgen import errors in the next task
        # (running gaussian on a different partition)
        if "prev_calc_molecule" in fw_spec:
            del fw_spec["prev_calc_molecule"]
        if "gaussian_input" in fw_spec:
            del fw_spec["gaussian_input"]
