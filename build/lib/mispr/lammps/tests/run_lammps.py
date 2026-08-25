import os

from fireworks import Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from mispr.lammps.fireworks.core import RunLammpsFW

ELEMENTS_LIST = ["N", "C", "C", "C", "H", "O", "O", "H", "O", "H", "Na"]
EMIN_SETTINGS = {
    "neigh_args": "2.0 bin",
    "neigh_modify_args": "delay 0 every 1 check yes",
    "pair_style": "lj/cut/coul/long",
    "pair_style_args": 10.0,
    "kspace_style": "pppm",
    "kspace_style_args": "1.0e-4",
    "data_file_name": "complex.data",
    "restart_final_name": "restart.emin.restart",
    "dump_file_name": "dump.emin.dump",
    "dump_modify_args": "",
}

NPT_SETTINGS = {
    "neigh_args": "2.0 bin",
    "neigh_modify_args": "delay 0 every 1 check yes",
    "special_bonds_style": "amber",
    "special_bonds_value": "",
    "bond_style": "harmonic",
    "bond_style_args": "",
    "angle_style": "harmonic",
    "dihedral_style": "harmonic",
    "improper_style": "cvff",
    "pair_style": "lj/cut/coul/long",
    "pair_style_args": 10.0,
    "kspace_style": "pppm",
    "kspace_style_args": "1.0e-4",
    "restart_filename": "restart.emin.restart",
    "group_definitions": "\n".join(
        ["group li 11", "group tfsi 1 2 3 4 5 6", "group solv 7 8 9 10"]
    ),
    "temperature_initial": 298.15,
    "velocity_seed": 250,
    "shake_logic": "#",
    "shake_group": "",
    "shake_topologies": "",
    "temperature_final": 298.15,
    "temp_damp": 100.0,
    "pressure_type": "iso",
    "pressure_initial": 1.0,
    "pressure_final": 1.0,
    "pres_damp": 1000.0,
    "thermo": 1000,
    "timestep": 1,
    "dump_period": 50000,
    "dump_file_name": "dump.npt.*.dump",
    "dump_modify_logic": "#",
    "dump_modify_args": "",
    "restart_period": 1000000,
    "restart_intermediate_file_name": "restart.npt.*.restart",
    "run": 2000000,
    "restart_final_file_name": "restart.npt.restart",
}

if __name__ == "__main__":
    # set up the LaunchPad and reset it
    launchpad = LaunchPad(
        host="mongodb+srv://mbliss01:idlewide@gettingstarted.dt0sv.mongodb.net/fireworks",
        uri_mode=True,
    )
    launchpad.reset("", require_password=False)

    # set working directory
    working_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "test_files", "run_lammps"
    )

    # settings
    control_filename = "test_npt.lammpsin"
    template_type = "npt"

    # define workflow
    spec = {}
    fw = RunLammpsFW(
        working_dir=working_dir,
        control_filename=control_filename,
        template_filename=template_type,
        control_settings=NPT_SETTINGS,
    )

    workflow = Workflow([fw])
    launchpad.add_wf(workflow)
    rapidfire(launchpad)
