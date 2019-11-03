# coding: utf-8


# This module parses Gaussian output files 

import os
import json

from fireworks import FiretaskBase, explicit_serialize, FWAction
from pymatgen import Molecule, IMolecule
from pymatgen.io.gaussian import GaussianOutput, GaussianInput

__author__ = 'Rasha Atwi'
__version__ = '0.1'
__email__ = 'rasha.atwi@tufts.edu'
__date__ = 'Aug 8, 2019'


@explicit_serialize
class ParseGaussianOutputFile(FiretaskBase):
    """
    This class reads a gaussian output file, converts it to an object, and saves its information to mongodb.
    """
    
    def run_task(self, fw_spec):
        files_dir = fw_spec["filesDir"]
        flag = 0
        outputs_dict = {}
        for gau_out in os.listdir(files_dir):
            if gau_out.endswith(".out") or gau_out.endswith(".log"):
                out = json.loads(json.dumps(GaussianOutput(f'{files_dir}/{gau_out}').as_dict()))
                del out['@module']
                del out['@class']
                name = os.path.splitext(gau_out)[0]
                outputs_dict[name] = out
                flag = 1
        if flag == 0:
            print("No Gaussian output files found")
        return FWAction(stored_data={'outputs': outputs_dict}, mod_spec=[{'_push': {'outputs': outputs_dict}}])

@explicit_serialize
class ConvertToXYZFile(FiretaskBase):
    """
    This class coverts a molecule object to xyz file.
    """
    def run_task(self, fw_spec):
        files_dir = fw_spec["filesDir"]
        var_dict = fw_spec[fw_spec['var']][0]
        temp_dict = {}
        if fw_spec['var'] == 'outputs':
            for key in var_dict:
                temp_dict[key] = var_dict[key]['output']['molecule']
            var_dict = temp_dict
        for name, var in var_dict.items():
            molecule = Molecule.from_sites(var)
            xyz_file = XYZ(molecule)
            xyz_file.write_file(f'{files_dir}/{name}.xyz')

# @explicit_serialize
# class ConvertGaussianOutputToInput(FiretaskBase):
#     """
#     This class converts a gaussian output to a gaussian input file.
#     """
#     def run_task(self, fw_spec):
#         outputs_dict = fw_spec['molecules'][0]
#         functional = fw_spec['functional']
#         basis_set = fw_spec['basis_set']
#         route_parameters = fw_spec['route_parameters']
#         link0_parameters = fw_spec['link0_parameters']
#         for name, out in outputs_dict.items():
#             molecule = Molecule.from_dict(out)
#             # TODO: automatic naming of the checkpoint files in link0_parameters
#             gau_in = molecule.to_input(molecule, title=name,
#                                    functional=functional, basis_set=basis_set, route_parameters=route_parameters,
#                                    link0_parameters=link0_parameters)
#             gau_in.write_file(f'/Users/rashaatwi/Desktop/{name}.com', cart_coords=True)
#             # TODO: find a more efficient way to delete last 4 blank lines of input file
#             # TODO: write an if function to either ask user for dir or get current dir
#             lines = open(f'/Users/rashaatwi/Desktop/{name}.com', 'r').readlines()
#             lines = lines[:-4]
#             open(f'/Users/rashaatwi/Desktop/{name}.com', 'w').writelines(lines)
#



