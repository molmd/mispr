import os

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from mispr.lammps.firetasks.run import RunTleap
from mispr.lammps.firetasks.write_inputs import WriteTleapScript

if __name__ == "__main__":

    # set up the LaunchPad and reset it
    launchpad = LaunchPad()
    launchpad.reset("", require_password=False)

    working_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "test_files", "tleap"
    )
    file_label = "dhps"

    tleap_settings = {
        "mol2_file_path": os.path.join(working_dir, f"{file_label}.mol2"),
        "frcmod_file_path": os.path.join(working_dir, f"{file_label}.frcmod"),
        "prmtop_file_path": os.path.join(working_dir, f"{file_label}.prmtop"),
        "inpcrd_file_path": os.path.join(working_dir, f"{file_label}.inpcrd"),
    }

    task1 = WriteTleapScript(working_dir=working_dir, tleap_settings=tleap_settings)

    task2 = RunTleap(working_dir=working_dir)

    # assemble FireWork from tasks and give the FireWork a unique id
    fire_work1 = Firework([task1, task2], name="RunTleap", fw_id=1)

    # assemble Workflow from FireWorks and their connections by id
    wf = Workflow([fire_work1])

    # store workflow and launch it
    launchpad.add_wf(wf)
    rapidfire(launchpad)
