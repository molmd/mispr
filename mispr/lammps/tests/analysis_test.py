import os
import sys

print(sys.path)

import infrastructure.lammps.fireworks.core_custom as ilfcc
from fireworks import LaunchPad
from fireworks.core.rocket_launcher import rapidfire

if __name__ == "__main__":
    print(sorted(sys.modules.keys()))
    # # set up the LaunchPad and reset it
    launchpad = LaunchPad(host="mongodb+srv://mbliss01:idlewide@gettingstarted.dt0sv.mongodb.net/fireworks",
                          uri_mode=True)
    launchpad.reset('', require_password=False)
