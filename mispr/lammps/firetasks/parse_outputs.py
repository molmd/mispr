# coding: utf-8

# Defines firetasks for reading output of Ambertools calcs and LAMMPS simulations

import os

from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from pymatgen.io.ambertools import PrmtopParser
from pymatgen.core.structure import Molecule

from mispr.lammps.utils.utils import add_ff_labels_to_dict
from mispr.gaussian.utilities.metadata import get_mol_formula

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@tufts.edu"
__status__ = "Development"
__date__ = "Apr 15, 2020"
__version__ = "0.0.1"


@explicit_serialize
class ProcessPrmtop(FiretaskBase):
    _fw_name = "Process Prmtop"
    required_params = ["molecule"]
    optional_params = [
        "working_dir",
        "prmtop_path",
        "prmtop_filename",
        "prmtop_dir",
        "unique_molecule_name",
        "system_force_field_dict",
    ]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        if isinstance(self.get("prmtop_path"), str):
            prmtop_file_path = self.get("prmtop_path")

        elif isinstance(self.get("prmtop_filename"), str):
            prmtop_filename = self.get("prmtop_filename")
            prmtop_dir = self.get("prmtop_dir", working_dir)
            prmtop_file_path = os.path.join(prmtop_dir, prmtop_filename)

        else:
            raise Exception(
                "Neither name nor path to a prmtop file was provided; "
                'either provide prmtop file name as "prmtop_filename" or '
                'path to prmtop as "prmtop_path".'
            )

        if not os.path.exists(prmtop_file_path):
            raise Exception(
                '"prmtop_file_path" is not a valid path; check that '
                '"prmtop_path" is correct or that "prmtop_filename" '
                'exists in "prmtop_dir"'
            )

        molecule = self.get("molecule")

        if not isinstance(molecule, Molecule):
            raise Exception('"molecule" is not a pymatgen Molecule object')

        # TODO: decide what unique id for the ff labels. If staying as smiles
        #  string, decide how to get smiles.
        unique_mol_name = self.get("unique_molecule_name", get_mol_formula(molecule))

        ff_param_dict_general = PrmtopParser(prmtop_file_path, molecule, "").to_dict()
        ff_param_dict_system = add_ff_labels_to_dict(
            ff_param_dict_general, unique_mol_name
        )

        # print(ff_param_dict_general)

        sys_ff_dict = fw_spec.get(
            "system_force_field_dict", self.get("system_force_field_dict", {})
        )
        sys_ff_dict[unique_mol_name] = ff_param_dict_system
        return FWAction(
            stored_data={
                "ff_param_dict_system": ff_param_dict_system,
                "ff_param_dict_general": ff_param_dict_general,
            },
            update_spec={"system_force_field_dict": sys_ff_dict},
        )
        # mod_spec = [{'_set': {"system_force_field_dict": sys_ff_dict}}])


if __name__ == "__main__":
    test_string = "test_string\t\n"
    print([test_string])
    print([test_string.strip()])
