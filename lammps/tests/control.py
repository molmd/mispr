import os
import pathlib
import numpy as np
from collections import OrderedDict

from fireworks import LaunchPad, Firework, FiretaskBase, FWAction, Workflow, FileWriteTask, FileTransferTask, explicit_serialize
from fireworks.core.rocket_launcher import rapidfire, launch_rocket
from infrastructure.lammps.fireworks.core_custom import AmbertoolsTasks
from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.ambertools import PrmtopParser
from infrastructure.lammps.firetasks.write_inputs import WriteDataFile, WriteControlFile, WriteTleapScript,\
    TLEAP_SETTING_DEFAULTS, TEMPLATE_DIR
from infrastructure.lammps.firetasks.run import RunLammps, RunAntechamber, RunParmchk, RunTleap
from infrastructure.lammps.firetasks.parse_outputs import ProcessPrmtop
from infrastructure.gaussian.utils.utils import get_mol_formula

ELEMENTS_LIST = ["N", "C", "C", "C", "H", "O", "O", "H", "O", "H", "Na"]

EMIN_SETTINGS = {"special_bonds_style": "amber",
                 "special_bonds_value": "",
                 "bond_style": "harmonic",
                 "bond_style_args": "",
                 "angle_style": "harmonic",
                 "dihedral_style": "harmonic",
                 "improper_style": "cvff",
                 "pair_style": "lj/cut/coul/long",
                 "pair_style_args": 10.0,
                 "kspace_style": "pppm",
                 "kspace_style_args": 1.0e-4,
                 "data_file_name": "complex.data",
                 "restart_final_name": "emin.restart",
                 "dump_modify_elements": ' '.join(ELEMENTS_LIST)}

if __name__ == "__main__":
    # set up the LaunchPad and reset it
    launchpad = LaunchPad(host="mongodb://superuser:idlewide@localhost:27017/fireworks?authSource=admin", uri_mode=True)
    launchpad.reset('', require_password=False)

    # set working directory
    working_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_files", "control")


    # # # Test Case 1: provide general template as string
    control_filename = "case_1.lammpsin"
    template_filename = "emin"
    template_path = os.path.join(TEMPLATE_DIR, template_filename)
    with open(template_path, 'r') as file:
        template_string = file.read()

    # # Define workflow for Case 1
    spec = {}
    t = [WriteControlFile(working_dir = working_dir,
                          control_filename = control_filename,
                          template_str = template_string,
                          control_settings = EMIN_SETTINGS)]

    firework_c1 = Firework(t,
                           spec = spec,
                           name = "WriteControlCase1",
                           fw_id = 1)

    wf_c1 = Workflow([firework_c1])

    # Store workflow for Case 1 and launch it
    launchpad.add_wf(wf_c1)
    rapidfire(launchpad)


    # # # Test Case 2: choose specific template file that in TEMPLATE_DIR
    launchpad.reset('', require_password=False)

    control_filename = "case_2.lammpsin"
    template_type = "emin"

    # # Define workflow for Case 2
    spec = {}
    t = [WriteControlFile(working_dir = working_dir,
                          control_filename = control_filename,
                          template_filename = template_type,
                          control_settings = EMIN_SETTINGS)]

    firework_c2 = Firework(t,
                           spec = spec,
                           name = "WriteControlCase2",
                           fw_id = 2)

    wf_c2 = Workflow([firework_c2])

    # Store workflow for Case 2 and launch it
    launchpad.add_wf(wf_c2)
    rapidfire(launchpad)


    # # # Test Case 3: provide general template as file
    launchpad.reset('', require_password=False)

    control_filename = "case_3.lammpsin"
    template_filename = "emin"
    template_path = TEMPLATE_DIR

    # # Define workflow for Case 3
    spec = {}
    t = [WriteControlFile(working_dir = working_dir,
                          control_filename = control_filename,
                          template_filename = template_filename,
                          template_dir = template_path,
                          control_settings = EMIN_SETTINGS)]

    firework_c3 = Firework(t,
                           spec = spec,
                           name = "WriteControlCase3",
                           fw_id = 3)

    wf_c3 = Workflow([firework_c3])

    # Store workflow for Case 3 and launch it
    launchpad.add_wf(wf_c3)
    rapidfire(launchpad)