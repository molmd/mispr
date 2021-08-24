import os

from collections import OrderedDict

import numpy as np

from fireworks import Firework, Workflow, LaunchPad, FiretaskBase, explicit_serialize
from fireworks.core.rocket_launcher import rapidfire

from pymatgen.io.gaussian import GaussianOutput
from pymatgen.core.structure import Molecule

from mispr.gaussian.utilities.metadata import get_mol_formula
from mispr.lammps.firetasks.parse_outputs import ProcessPrmtop


@explicit_serialize
class PrintFW(FiretaskBase):
    """
    Firetask for confirming that modspec works as intended in ProcessPrmtop firetask
    """

    def run_task(self, fw_spec):
        print(str(fw_spec["system_force_field_dict"]))


if __name__ == "__main__":

    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset("", require_password=False)

    working_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "test_files", "prmtop"
    )

    # write inputs for initial
    # Phen_type parameter prep
    Phen_type = "dhps"
    Phen_type_gaussian_output = GaussianOutput(
        os.path.join(working_dir, Phen_type + ".out")
    )
    Phen_type_molecule = Phen_type_gaussian_output.structures[-1]
    Phen_type_molecule.set_charge_and_spin(
        Phen_type_gaussian_output.charge, Phen_type_gaussian_output.spin_multiplicity
    )
    Phen_type_label = get_mol_formula(Phen_type_molecule)
    # Spce parameter prep
    Spce_molecule = Molecule.from_file(os.path.join(working_dir, "SPC_E.pdb"))
    Spce_label = get_mol_formula(Spce_molecule)
    Spce_param_dict = {
        "Molecule": Spce_molecule,
        "Labels": ["ow" + Spce_label, "hw" + Spce_label, "hw" + Spce_label],
        "Masses": OrderedDict({"ow" + Spce_label: 16.000, "hw" + Spce_label: 1.008}),
        "Nonbond": [[0.155394259, 3.16555789], [0.0, 0.0]],
        "Bonds": [
            {"coeffs": [553.0, 1], "types": [("ow" + Spce_label, "hw" + Spce_label)]}
        ],
        "Angles": [
            {
                "coeffs": [100.0, 109.47],
                "types": [("hw" + Spce_label, "ow" + Spce_label, "hw" + Spce_label)],
            }
        ],
        "Dihedrals": [],
        "Impropers": [],
        "Improper Topologies": None,
        "Charges": np.asarray([-0.8476, 0.4238, 0.4238]),
    }
    # Na parameter prep
    Na_molecule = Molecule.from_file(os.path.join(working_dir, "Na.pdb"))
    Na_molecule.set_charge_and_spin(1)
    Na_label = get_mol_formula(Na_molecule)
    Na_param_dict = {
        "Molecule": Na_molecule,
        "Labels": ["na" + Na_label],
        "Masses": OrderedDict({"na" + Na_label: 22.99}),
        "Nonbond": [[0.02639002, 2.590733472]],  # from frcmod.ions1lm_126_spce (2015)
        #                     'Nonbond': [[0.3526418, 2.159538]], # from frcmod.ionsjc_tip3p (2008)
        "Bonds": [],
        "Angles": [],
        "Dihedrals": [],
        "Impropers": [],
        "Improper Topologies": None,
        "Charges": np.asarray([1.0]),
    }

    sys_ff_dict = {Na_label: Na_param_dict, Spce_label: Spce_param_dict}

    spec = {"system_force_field_dict": sys_ff_dict}

    task1 = ProcessPrmtop(
        working_dir=working_dir,
        molecule=Phen_type_molecule,
        prmtop_filename="dhps.prmtop",
    )

    task2 = PrintFW()

    # assemble FireWork from tasks and give the FireWork a unique id
    fire_work1 = Firework([task1, task2], spec=spec, name="ParsePrmtop", fw_id=1)

    # assemble Workflow from FireWorks and their connections by id
    wf = Workflow([fire_work1])

    # store workflow and launch it
    launchpad.add_wf(wf)
    rapidfire(launchpad)
