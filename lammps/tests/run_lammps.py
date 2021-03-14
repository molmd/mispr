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
from infrastructure.lammps.fireworks.core_custom import RunLammpsFW

ELEMENTS_LIST = ["N", "C", "C", "C", "H", "O", "O", "H", "O", "H", "Na"]
EMIN_SETTINGS = {"neigh_args": "2.0 bin",
                 "neigh_modify_args": "delay 0 every 1 check yes",
                 "pair_style": "lj/cut/coul/long",
                 "pair_style_args": 10.0,
                 "kspace_style": "pppm",
                 "kspace_style_args": "1.0e-4",
                 "data_file_name": "complex.data",
                 "restart_final_name": "restart.emin.restart",
                 "dump_file_name": "dump.emin.dump",
                 "dump_modify_args": ""}

if __name__ == "__main__":
    # set up the LaunchPad and reset it
    launchpad = LaunchPad(host="mongodb://superuser:idlewide@localhost:27017/fireworks?authSource=admin", uri_mode=True)
    launchpad.reset('', require_password=False)

    # set working directory
    working_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_files", "run_lammps")

    # settings
    control_filename = "test.lammpsin"
    template_type = "emin_gaff"

    # # Define workflow
    spec = {}
    fw = RunLammpsFW(working_dir = working_dir,
                     control_filename = control_filename,
                     template_filename = template_type,
                     control_settings = EMIN_SETTINGS)

    workflow = Workflow([fw])
    launchpad.add_wf(workflow)
    rapidfire(launchpad)