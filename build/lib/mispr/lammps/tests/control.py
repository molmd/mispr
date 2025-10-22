import os

from fireworks import Firework, Workflow, LaunchPad
from fireworks.core.rocket_launcher import rapidfire

from mispr.lammps.firetasks.write_inputs import WriteControlFile

ELEMENTS_LIST = ["N", "C", "C", "C", "H", "O", "O", "H", "O", "H", "Na"]

EMIN_SETTINGS = {
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
    "kspace_style_args": 1.0e-4,
    "data_file_name": "complex.data",
    "restart_final_name": "emin.restart",
    "dump_modify_elements": " ".join(ELEMENTS_LIST),
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
    "group_definitions": [
        "group li 11",
        "group tfsi 1 2 3 4 5 6",
        "group solv 7 8 9 10",
    ],
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
        os.path.dirname(os.path.abspath(__file__)), "test_files", "control"
    )

    # Test Case 1: provide general template as string
    # control_filename = "case_1.lammpsin"
    # template_filename = "emin"
    # template_path = os.path.join(TEMPLATE_DIR, template_filename)
    # with open(template_path, 'r') as file:
    #     template_string = file.read()
    #
    # # # Define workflow for Case 1
    # spec = {}
    # t = [WriteControlFile(working_dir = working_dir,
    #                       control_filename = control_filename,
    #                       template_str = template_string,
    #                       control_settings = EMIN_SETTINGS)]
    #
    # firework_c1 = Firework(t,
    #                        spec = spec,
    #                        name = "WriteControlCase1",
    #                        fw_id = 1)
    # wf_c1 = Workflow([firework_c1])
    # Store workflow for Case 1 and launch it
    # launchpad.add_wf(wf_c1)
    # rapidfire(launchpad)

    # Test Case 2: choose specific template file that in TEMPLATE_DIR
    launchpad.reset("", require_password=False)

    control_filename = "case_2_npt_group_list.lammpsin"
    template_type = "npt"

    # # Define workflow for Case 2
    spec = {}
    t = [
        WriteControlFile(
            working_dir=working_dir,
            control_filename=control_filename,
            template_filename=template_type,
            control_settings=NPT_SETTINGS,
        )
    ]

    firework_c2 = Firework(t, spec=spec, name="WriteControlCase2", fw_id=2)
    wf_c2 = Workflow([firework_c2])
    # Store workflow for Case 2 and launch it
    launchpad.add_wf(wf_c2)
    rapidfire(launchpad)

    # Test Case 3: provide general template as file
    # launchpad.reset('', require_password=False)
    #
    # control_filename = "case_3.lammpsin"
    # template_filename = "emin"
    # template_path = TEMPLATE_DIR
    #
    # # # Define workflow for Case 3
    # spec = {}
    # t = [WriteControlFile(working_dir = working_dir,
    #                       control_filename = control_filename,
    #                       template_filename = template_filename,
    #                       template_dir = template_path,
    #                       control_settings = EMIN_SETTINGS)]
    #
    # firework_c3 = Firework(t,
    #                        spec = spec,
    #                        name = "WriteControlCase3",
    #                        fw_id = 3)
    #
    # wf_c3 = Workflow([firework_c3])
    # Store workflow for Case 3 and launch it
    # launchpad.add_wf(wf_c3)
    # rapidfire(launchpad)
