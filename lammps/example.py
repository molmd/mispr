from fireworks import LaunchPad, Firework, FiretaskBase, FWAction, Workflow
from fireworks.core.rocket_launcher import rapidfire, launch_rocket
from infrastructure.lammps.fireworks.core_custom import AmbertoolsTasks

if __name__ == "__main__":

    # set up the LaunchPad and reset it
    launchpad = LaunchPad()