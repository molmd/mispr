import os

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from mispr.lammps.firetasks.run import RunParmchk

if __name__ == "__main__":

    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset("", require_password=False)

    working_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "test_files", "parmchk"
    )

    input_filename = "dhps.mol2"
    output_filename = "dhps.frcmod"

    task1 = RunParmchk(
        working_dir=working_dir,
        input_filename_p=input_filename,
        output_filename_p=output_filename,
    )

    # assemble FireWork from tasks and give the FireWork a unique id
    fire_work1 = Firework(task1, name="RunParm", fw_id=1)

    # assemble Workflow from FireWorks and their connections by id
    wf = Workflow([fire_work1])

    # store workflow and launch it
    launchpad.add_wf(wf)
    rapidfire(launchpad)
