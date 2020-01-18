print('Write Inputs')
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
__version__ = '0.1'
__email__ = 'rasha.atwi@tufts.edu'
__date__ = 'Aug 8, 2019'


@explicit_serialize
class WriteInput(FiretaskBase):
    required_params = []
    optional_params = ["gaussian_input", "molecule", "gaussian_input_params",
                       "input_file", "write_to_dir", "cart_coords",
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
                 'charge':  mol_copy.charge}

    def run_task(self, fw_spec):
        # TODO: PASS A MOLECULE FROM THE DATABASE
        input_file = os.path.join(self.get("write_to_dir", ""),
                                  self.get("input_file", "mol.com"))

        # if a full Gaussian object is provided
        if self.get("gaussian_input") and isinstance(self.get("gaussian_input"),
                                                     GaussianInput):
            gaussin = self["gaussian_input"]
        # if a molecule is being passed through fw_spec
        elif fw_spec.get("prev_calc_molecule"):
            prev_calc_mol = fw_spec.get("prev_calc_molecule")
            # if a molecule is also passed as an optional parameter
            if self.get("molecule"):
                mol = self.get("molecule")
                # check if mol and prev_calc_mol are isomorphic
                mol_graph = MoleculeGraph.with_local_env_strategy(mol,
                                                                  OpenBabelNN(),
                                                                  reorder=False,
                                                                  extend_structure=False)
                prev_mol_graph = MoleculeGraph.with_local_env_strategy(prev_calc_mol,
                                                                       OpenBabelNN(),
                                                                       reorder=False,
                                                                       extend_structure=False)
                # If they are isomorphic, aka a previous FW has not changed bonding,
                # then we will use prev_calc_mol. If bonding has changed, we will use mol.
                if mol_graph.isomorphic_to(prev_mol_graph):
                    mol = prev_calc_mol
                else:
                    print("Not using prev_calc_mol as it is not isomorphic to passed molecule!")
            else:
                mol = prev_calc_mol
            self._update_charge(mol)
            gaussin = GaussianInput(mol, **self.get("gaussian_input_params", {}))
        # if a molecule is only included as an optional parameter
        elif self.get("molecule"):
            self._update_charge(self.get("molecule"))
            gaussin = GaussianInput(self.get("molecule"), **self.get("gaussian_input_params", {}))
        # if no molecule is present raise an error
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )
        gaussin.write_file(input_file, self.get("cart_coords", False))
        fw_spec['input_file_path'] = input_file
        input_dict = gaussin.as_dict()
        # fw_spec['run'].update({'input': input_dict })
        del fw_spec['prev_calc_molecule']
        print(fw_spec)
        return FWAction(update_spec=fw_spec)


@explicit_serialize
class ConvertMoleculeToGaussianInput(FiretaskBase):
    # TODO: allow taking a mol object from database directly and/or using the above class
    def run_task(self, fw_spec):
        files_dir = fw_spec["filesDir"]
        var_dict = fw_spec[fw_spec['var']][0]
        temp_dict = {}
        if fw_spec['var'] == 'outputs':
            for key in var_dict:
                temp_dict[key] = var_dict[key]['output']['molecule']
            var_dict = temp_dict
        functional = fw_spec['functional']
        basis_set = fw_spec['basis_set']
        route_parameters = fw_spec['route_parameters']
        link0_parameters = fw_spec['link0_parameters']
        input_parameters = fw_spec['input_parameters']
        # TODO: define job names for hybrid jobs
        # TODO: remove job names from previous calculation
        if 'Opt' in route_parameters:
            job = 'Opt'
        elif 'Freq' in route_parameters:
            job = 'Freq'
        elif 'NMR' in route_parameters:
            job = 'NMR'
        elif 'pop' in route_parameters and 'ESP' in input_parameters:
            job = 'ESP'
        for name, var in var_dict.items():
            molecule = Molecule.from_sites(var)
            molecule_copy = deepcopy(molecule)
            # TODO: find a way to calculate charge correctly
            molecule_copy.add_oxidation_state_by_element(
                {"Mg": 2, "Cl": -1, "N": -1, "S": 0, "O": 0, "F": 0, "C": 0, "H": 0})
            molecule_copy.set_charge_and_spin(super(IMolecule, molecule_copy).charge)
            # TODO: automatic naming of the checkpoint files in link0_parameters
            gau_in = GaussianInput(molecule, charge=molecule_copy.charge,
                                   spin_multiplicity=None, title=None,
                                   functional=functional, basis_set=basis_set,
                                   route_parameters=route_parameters,
                                   link0_parameters=link0_parameters,
                                   input_parameters=input_parameters)
            gau_in.write_file(f'{files_dir}/{name}_{job}.com', cart_coords=True)
            # TODO: find a more efficient way to delete last 4 blank lines of input file
            if not input_parameters:
                lines = open(f'{files_dir}/{name}_{job}.com', 'r').readlines()
                lines = lines[:-2]
                open(f'{files_dir}/{name}_{job}.com', 'w').writelines(lines)
