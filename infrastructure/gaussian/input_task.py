# coding: utf-8


# This module defines firetasks for writing Gaussian input files

import os
import subprocess
from copy import deepcopy

from fireworks import FiretaskBase, explicit_serialize, FWAction
from fireworks.core.rocket_launcher import rapidfire
from pymatgen import Molecule, IMolecule
from pymatgen.io.gaussian import GaussianInput

__author__ = 'Rasha Atwi'
__version__ = '0.1'
__email__ = 'rasha.atwi@tufts.edu'
__date__ = 'Aug 8, 2019'


# Omitting _fw_name and identifying the class by the module name and class name
@explicit_serialize
class ConvertToMoleculesTask(FiretaskBase):
    """
    This class reads a molecule from a file, converts it to a mol object, and saves it as dict to mongodb.
    Supported formats include
        xyz|pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv,
        gaussian input (gjf|g03|g09|com|inp),
        Gaussian output (.out), and
        pymatgen's JSON serialized molecules.
    Requires openbabel to be installed.
    """
    # _fw_name = "Convert To Molecules Task"

    def run_task(self, fw_spec):
        files_dir = fw_spec["filesDir"]
        mols_dict = {}
        flag = 0
        for file in os.listdir(files_dir):
            mol = Molecule.from_file(f'{files_dir}/{file}').as_dict()
            name = os.path.splitext(file)[0]
            mols_dict[name] = mol
            flag = 1
        if flag == 0:
            print("No molecule files found")
        return FWAction(stored_data={'molecules': mols_dict}, mod_spec=[{'_push': {'molecules': mols_dict}}])


@explicit_serialize
class ConvertMoleculeToGaussianInputTask(FiretaskBase):
    """
    This class inherits FireTaskBase class; all methods of FireTaskBase become methods of this class.

    """
    # _fw_name = "Convert Molecule To Gaussian Input Task"

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
        #elif 'pop' in route_parameters and f'{name}.esp' in input_parameters:
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
            gau_in = GaussianInput(molecule, charge=molecule_copy.charge, spin_multiplicity=None, title=None,
                                   functional=functional, basis_set=basis_set, route_parameters=route_parameters,
                                   link0_parameters=link0_parameters, input_parameters=input_parameters)
            gau_in.write_file(f'{files_dir}/{name}_{job}.com', cart_coords=True)
            # TODO: find a more efficient way to delete last 4 blank lines of input file
            if not input_parameters:
                lines = open(f'{files_dir}/{name}_{job}.com', 'r').readlines()
                lines = lines[:-2]
                open(f'{files_dir}/{name}_{job}.com', 'w').writelines(lines)

@explicit_serialize
class RunGaussianDirect(FiretaskBase):
    """
    This class executes a command directly (no custodian)
    """
    # TODO: replace with a custodian
    def run_task(self, fw_spec):
        files_dir = fw_spec["filesDir"]
        flag = 0
        for in_file in os.listdir(files_dir):
            if in_file.endswith(".com") or in_file.endswith(".in"):
                name = os.path.splitext(in_file)[0]
                subprocess.call(['g09', name + '.com', name + '.out'])
                flag = 1
        if flag == 0:
            print("No Gaussian input file found")

