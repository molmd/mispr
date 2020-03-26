# coding: utf-8


# Defines firetasks for writing Gaussian input files.

import os
from copy import deepcopy

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.core.structure import Molecule, IMolecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.io.gaussian import GaussianInput

__author__ = 'Rasha Atwi'
__version__ = '0.0.1'
__email__ = 'rasha.atwi@tufts.edu'
__date__ = 'Aug 8, 2019'


@explicit_serialize
class WriteInput(FiretaskBase):
    _fw_name = "Write Gaussian Input File"

    required_params = []
    optional_params = ["gaussian_input", "molecule", "gaussian_input_params",
                       "input_file", "working_dir", "cart_coords",
                       "oxidation_states"]

    def _update_charge(self, mol):
        """
        Calculates the charge of a molecule/cluster using the oxidation state
        of individual elements.
        """
        if self.get("oxidation_states"):
            mol_copy = deepcopy(mol)
            mol_copy.add_oxidation_state_by_element(
                self.get("oxidation_states"))
            mol_copy.set_charge_and_spin(super(IMolecule, mol_copy).charge)
            self["gaussian_input_params"] = \
                {**self.get("gaussian_input_params", {}),
                 'charge': mol_copy.charge}

    def run_task(self, fw_spec):
        working_dir = fw_spec. \
            get('working_dir', self.get('working_dir', os.getcwd()))

        input_file = self.get("input_file", "mol.com")
        input_path = os.path.join(working_dir, input_file)

        # if a full Gaussian object is provided
        if self.get("gaussian_input") and \
                isinstance(self.get("gaussian_input"), GaussianInput):
            gaussin = self["gaussian_input"]

        # if a full Gaussian object is being passed through fw_spec
        elif fw_spec.get("gaussian_input") and \
                isinstance(fw_spec.get("gaussian_input"), GaussianInput):
            gaussin = fw_spec.get("gaussian_input")

        # if a molecule is being passed through fw_spec
        elif fw_spec.get("prev_calc_molecule"):
            prev_calc_mol = fw_spec.get("prev_calc_molecule")
            # if a molecule is also passed as an optional parameter
            if self.get("molecule"):
                mol = self.get("molecule")
                # check if mol and prev_calc_mol are isomorphic
                mol_graph = MoleculeGraph. \
                    with_local_env_strategy(mol, OpenBabelNN(), reorder=False,
                                            extend_structure=False)
                prev_mol_graph = MoleculeGraph. \
                    with_local_env_strategy(prev_calc_mol, OpenBabelNN(),
                                            reorder=False,
                                            extend_structure=False)
                # If they are isomorphic, aka a previous FW has not changed bonding,
                # then we will use prev_calc_mol. If bonding has changed, we will use mol.
                if mol_graph.isomorphic_to(prev_mol_graph):
                    mol = prev_calc_mol
                else:
                    print(
                        "Not using prev_calc_mol as it is not isomorphic to passed molecule!")
            else:
                mol = prev_calc_mol
            self._update_charge(mol)
            gaussin = \
                GaussianInput(mol, **self.get("gaussian_input_params", {}))
        # if a molecule is only included as an optional parameter
        elif self.get("molecule"):
            self._update_charge(self.get("molecule"))
            gaussin = GaussianInput(self.get("molecule"),
                                    **self.get("gaussian_input_params", {}))
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

# TODO: write a firetask for creating an input file using default parameters
# and another one using custom ones that offer more flexibility
