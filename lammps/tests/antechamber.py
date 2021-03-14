import os

from fireworks import LaunchPad, Firework, FiretaskBase, FWAction, Workflow
from fireworks.core.rocket_launcher import rapidfire, launch_rocket
from infrastructure.lammps.fireworks.core_custom import AmbertoolsTasks
from pymatgen.core.structure import Molecule
from infrastructure.lammps.firetasks.write_inputs import WriteDataFile, WriteControlFile, WriteTleapScript,\
    TLEAP_SETTING_DEFAULTS
from infrastructure.lammps.firetasks.run import RunLammps, RunAntechamber, RunParmchk, RunTleap
from infrastructure.lammps.firetasks.parse_outputs import ProcessPrmtop
from infrastructure.gaussian.utils.utils import get_mol_formula

if __name__ == "__main__":

    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset('', require_password=False)

    working_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_files", "antechamber")
    # print(working_dir)
    input_filename = "dhps.esp"
    output_filename = "dhps.mol2"

    task1 = RunAntechamber(working_dir = working_dir,
                           input_filename_a = input_filename,
                           output_filename_a = output_filename)

    # assemble FireWork from tasks and give the FireWork a unique id
    fire_work1 = Firework(task1, name='RunAnt', fw_id=1)

    # assemble Workflow from FireWorks and their connections by id
    wf = Workflow([fire_work1])

    # store workflow and launch it
    launchpad.add_wf(wf)
    rapidfire(launchpad)