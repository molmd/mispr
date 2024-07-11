import os

# default rdf settings
RDF_SETTINGS = {
    "rdf_type": "atomic",
    "r_cut": 20,
    "bin_size": [0.05],
    "filename": os.path.abspath(
        os.path.join("../../../lammps", "nvt", "dump.nvt.*.dump")
    ),
    "path_or_buff": "rdf.csv",
    "save_mode": True,
}

# default cn settings
CN_SETTINGS = {
    "cn_type": "atomic",
    "bin_size": [0.05],
    "filename": os.path.abspath(
        os.path.join("../../../lammps", "nvt", "dump.nvt.*.dump")
    ),
    "path_or_buff": "cn.csv",
    "save_mode": True,
}

# default cluster analysis settings
CLUSTERS_SETTINGS = {
    "full_trajectory": True,
    "alter_atom_types": False,
    "max_force": 0.75,
    "filename": os.path.abspath(
        os.path.join("../../../lammps", "nvt", "dump.nvt.*.dump")
    ),
    "find_top": True,
    "perc": None,
    "cum_perc": 90,
    "zip": True,
}

# default msd settings
MSD_SETTINGS = {
    "timestep": 1,
    "units": "real",
    "tao_coeff": 4,
    "msd_type": "com",
    "return_all": False,
    "com_drift": False,
    "avg_interval": False,
    "save_msd": True,
    "msd_method": "from_dump",
}

# default diffusion settings
DIFF_SETTINGS = {
    "initial_time": None,
    "final_time": None,
    "dimension": 3,
    "save": True,
    "plot": True,
    "diff_dist": False,
}

# force field dictionary
FF_DICT_KEYS = [
    "Molecule",
    "Labels",
    "Masses",
    "Nonbond",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
    "Improper Topologies",
    "Charges",
]

# lammps recipt template files
# TODO: decide if I want specific templates for melt and quench
TEMPLATE_TYPES = [
    "emin_general",
    "emin_gaff",
    "npt",
    "nvt",
    "min_interface",
    "nve_interface",
    "nvt_interface",
]

# default tleap settings
TLEAP_SETTINGS = {
    "source_file_path": "leaprc.gaff",
    "mol2_file_path": "mol.mol2",
    "frcmod_file_path": "mol.frcmod",
    "prmtop_file_path": "mol.prmtop",
    "inpcrd_file_path": "mol.inpcrd",
}

# default emin settings
EMIN_SETTINGS = {
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
    "pair_style_args": 12.0,
    "pair_modify_key": "mix",
    "pair_modify_value": "arithmetic",
    "kspace_style": "pppm",
    "kspace_style_args": "1.0e-4",
    "data_filename": "data.mixture",
    "restart_finalname": "restart.emin",
    "dump_filename": "dump.emin",
    "dump_modify_args": "",
    "data_final_filename": "data.emin",
}

# default npt settings
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
    "pair_style_args": 12.0,
    "pair_modify_key": "mix",
    "pair_modify_value": "arithmetic",
    "kspace_style": "pppm",
    "kspace_style_args": "1.0e-4",
    "restart_filename": "restart.emin",
    "group_definitions": "",
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
    "dump_filename": "dump.npt.*.dump",
    "dump_modify_logic": "#",
    "dump_modify_args": "",
    "restart_period": 1000000,
    "restart_intermediate_filename": "restart.npt_r",
    "run": 2000000,
    "restart_final_filename": "restart.npt",
    "data_final_filename": "data.npt",
}

# default nvt settings
NVT_SETTINGS = {
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
    "pair_style_args": 12.0,
    "pair_modify_key": "mix",
    "pair_modify_value": "arithmetic",
    "kspace_style": "pppm",
    "kspace_style_args": "1.0e-4",
    "restart_filename": "restart.npt",
    "group_definitions": "",
    "shake_logic": "#",
    "shake_group": "",
    "shake_topologies": "",
    "velocity_logic": "#",
    "velocity_seed": 250,
    "temperature_initial": 298.15,
    "temperature_final": 298.15,
    "temp_damp": 100.0,
    "thermo": 5,
    "compute_definitions": "",
    "thermo_style_compute": "",
    "timestep": 1,
    "reset_timestep_logic": "",
    "reset_timestep_value": 0,
    "dump_period": 50000,
    "dump_filename": "dump.nvt.*.dump",
    "dump_modify_logic": "#",
    "dump_modify_args": "",
    "restart_period": 1000000,
    "restart_intermediate_filename": "restart.nvt_r",
    "run": 5000000,
    "restart_final_filename": "restart.nvt",
    "data_final_filename": "data.nvt",
}

# default nve settings for the nve_interface template
MIN_INTERFACE_SETTINGS = {
    "log_filename": "lammps.log",
    "units_type": "real",
    "x_pbc": "p",
    "y_pbc": "p",
    "z_pbc": "f",
    "atom_style": "full",
    "n_dimension": 3,
    "neigh_args": "2.0 bin",
    "neigh_modify_args": "delay 0 every 1 check yes",
    "bond_style": "harmonic",
    "bond_style_values": "",
    "angle_style": "harmonic",
    "dihedral_style": "harmonic",
    "improper_style": "cvff",
    "pair_style": "lj/cut/coul/long",
    "pair_style_values": 12.0,
    "pair_modify_args": "mix arithmetic shift no tail no",
    "kspace_style": "pppm",
    "kspace_style_values": "1.0e-4",
    "kspace_modify_line": "kspace_modify slab 3.0",
    "data_filename": "interface.data",
    "atom_defniitions": "#no additional atom definitions used",
    "group_definitions_additional": "#no additional group definitions used",
    "run_style_args": "verlet",
    "thermo": 1,
    "thermo_style_additional": "#no additional thermo_style values used",
    "thermo_modify_args": "flush yes lost warn",
    "dump_logic": "#",
    "dump_period": 1,
    "dump_filename": "dump.min_interface.*.dump",
    "dump_attributes": "id mol type element mass q x y z xu yu zu ix iy iz vx vy vz fx fy fz",
    "dump_modify_logic": "#",
    "dump_modify_args": "",
    "dump_additional": "#no additional dump commands used",
    "setforce_logic": "",
    "setforce_args": "0.0 0.0 0.0",
    "fix_definitions_additional": "#no additional fix definitions used",
    "rattle_logic": "#",
    "rattle_tol": "1.0e-4",
    "rattle_iter": 200,
    "rattle_N": 0,
    "shake_logic": "#",
    "shake_tol": "1.0e-4",
    "shake_iter": 200,
    "shake_N": 0,
    "min_style_1": "sd",
    "minimize_args_1": "1e-6 1000 10 10000",
    "min_logic_2": "",
    "min_style_2": "cg",
    "minimize_args_2": "1e-6 100 99990 990000",
    "unfix_additional": "#no additional unfix commands used",
    "restart_final_filename": "restart.min_interface.restart",
    "data_final_filename": "data.min_interface.data",
    "data_final_args": "pair ij",
}

NVE_INTERFACE_SETTINGS = {
    "log_filename": "lammps.log",
    "units_type": "real",
    "x_pbc": "p",
    "y_pbc": "p",
    "z_pbc": "f",
    "atom_style": "full",
    "n_dimension": 3,
    "neigh_args": "2.0 bin",
    "neigh_modify_args": "delay 0 every 1 check yes",
    "bond_style": "harmonic",
    "bond_style_values": "",
    "angle_style": "harmonic",
    "dihedral_style": "harmonic",
    "improper_style": "cvff",
    "pair_style": "lj/cut/coul/long",
    "pair_style_values": 12.0,
    "pair_modify_args": "mix arithmetic shift no tail no",
    "kspace_style": "pppm",
    "kspace_style_values": "1.0e-4",
    "kspace_modify_line": "kspace_modify slab 3.0",
    "read_data_logic": "#",
    "data_filename": "interface.data",
    "read_restart_logic": "",
    "restart_filename": "restart.min_interface.restart",
    "atom_defniitions": "#no additional atom definitions used",
    "group_definitions_additional": "#no additional group definitions used",
    "timestep": 1,
    "run_style_args": "verlet",
    "compute_definitions_additional": "#no additional compute definitions used",
    "thermo": 1000,
    "thermo_style_additional": "#no additional thermo_style values used",
    "thermo_modify_args": "flush yes lost warn",
    "dump_period": 50000,
    "dump_filename": "dump.nve_interface.*.dump",
    "dump_attributes": "id mol type element mass q x y z xu yu zu ix iy iz vx vy vz fx fy fz",
    "dump_modify_logic": "#",
    "dump_modify_args": "",
    "dump_additional": "#no additional dump commands used",
    "restart_period": 1000000,
    "restart_intermediate_filename": "restart.nve_interface.*.restart",
    "temperature_initial": 298.15,
    "velocity_seed": 250,
    "velocity_args": "mom yes rot yes dist gaussian",
    "velocity_definitions_additional": "#no additional velocity definitions used",
    "setforce_logic": "",
    "setforce_args": "0.0 0.0 0.0",
    "rescale_steps": 10,
    "temperature_final": 298.15,
    "rescale_window": 5.0,
    "rescale_fraction": 1.0,
    "fix_rescale_modify_args": "Tfluid",
    "fix_definitions_additional": "#no additional fix definitions used",
    "rattle_logic": "#",
    "rattle_tol": "1.0e-4",
    "rattle_iter": 200,
    "rattle_N": 0,
    "shake_logic": "#",
    "shake_tol": "1.0e-4",
    "shake_iter": 200,
    "shake_N": 0,
    "reset_timestep_logic": "",
    "reset_timestep_args": 0,
    "run": 100000,
    "unfix_additional": "#no additional unfix commands used",
    "restart_final_filename": "restart.nve_interface.restart",
    "data_final_filename": "data.nve_interface.data",
    "data_final_args": "pair ij",
}

# default nve settings for the nvt_interface template
NVT_INTERFACE_SETTINGS = {
    "log_filename": "lammps.log",
    "units_type": "real",
    "x_pbc": "p",
    "y_pbc": "p",
    "z_pbc": "f",
    "atom_style": "full",
    "n_dimension": 3,
    "neigh_args": "2.0 bin",
    "neigh_modify_args": "delay 0 every 1 check yes",
    "bond_style": "harmonic",
    "bond_style_values": "",
    "angle_style": "harmonic",
    "dihedral_style": "harmonic",
    "improper_style": "cvff",
    "pair_style": "lj/cut/coul/long",
    "pair_style_values": 12.0,
    "pair_modify_args": "mix arithmetic shift no tail no",
    "kspace_style": "pppm",
    "kspace_style_values": "1.0e-4",
    "kspace_modify_line": "kspace_modify slab 3.0",
    "restart_filename": "restart.nve_interface.restart",
    "atom_defniitions": "#no additional atom definitions used",
    "group_definitions_additional": "#no additional group definitions used",
    "timestep": 1,
    "run_style_args": "verlet",
    "compute_definitions_additional": "#no additional compute definitions used",
    "thermo": 5,
    "thermo_style_additional": "#no additional thermo_style values used",
    "thermo_modify_args": "flush yes lost warn",
    "dump_period": 50000,
    "dump_filename": "dump.nvt_interface.*.dump",
    "dump_attributes": "id mol type element mass q x y z xu yu zu ix iy iz vx vy vz fx fy fz",
    "dump_modify_logic": "#",
    "dump_modify_args": "",
    "dump_additional": "#no additional dump commands used",
    "restart_period": 1000000,
    "restart_intermediate_filename": "restart.nvt_interface.*.restart",
    "temperature_initial": 298.15,
    "velocity_logic": "#",
    "velocity_seed": 250,
    "velocity_args": "mom yes rot yes dist gaussian",
    "velocity_definitions_additional": "#no additional velocity definitions used",
    "setforce_logic": "",
    "setforce_args": "0.0 0.0 0.0",
    "temperature_final": 298.15,
    "temp_damp": 100.0,
    "nvt_args_additional": "#no additional nvt arguments used",
    "fix_nvt_modify_args": "Tfluid",
    "fix_definitions_additional": "#no additional fix definitions used",
    "rattle_logic": "#",
    "rattle_tol": "1.0e-4",
    "rattle_iter": 200,
    "rattle_N": 0,
    "shake_logic": "#",
    "shake_tol": "1.0e-4",
    "shake_iter": 200,
    "shake_N": 0,
    "reset_timestep_logic": "#",
    "reset_timestep_args": "0",
    "run": 10000000,
    "unfix_additional": "#no additional unfix commands used",
    "restart_final_filename": "restart.nvt_interface.restart",
    "data_final_filename": "data.nvt_interface.data",
    "data_final_args": "pair ij",
}

# default lammps recipe
LAMMPS_RECIPE = [
    ["emin", ["template_filename", "emin_gaff"]],
    ["npt", ["template_filename", "npt"]],
    ["melt", ["template_filename", "nvt"]],
    ["quench", ["template_filename", "nvt"]],
    ["nvt", ["template_filename", "nvt"]],
]

# default lammps recipe settings
RECIPE_SETTINGS = [
    {"data_filename": "../data.mixture", "restart_finalname": "restart.emin"},
    {
        "restart_filename": "../emin/restart.emin",
        "restart_final_filename": "restart.npt",
    },
    {
        "restart_filename": "../npt/restart.npt",
        "temperature_initial": 500.0,
        "temperature_final": 500.0,
        "run": 2000000,
        "restart_final_filename": "restart.melt_500K",
        "data_final_filename": "data.melt",
    },
    {
        "restart_filename": "../melt/restart.melt_500K",
        "temperature_initial": 500.0,
        "temperature_final": 298.15,
        "run": 3000000,
        "restart_final_filename": "restart.quench_298K",
        "data_final_filename": "data.quench_298K",
    },
    {
        "restart_filename": "../quench/restart.quench_298K",
        "temperature_initial": 298.15,
        "temperature_final": 298.15,
        "run": 5000000,
        "restart_final_filename": "restart.nvt_5ns",
        "data_final_filename": "data.nvt_5ns",
    },
]

# default lammps run resources
QADAPTER_RUN_LAMMPS_SPEC = [
    {"walltime": "00:10:00", "job_name": "emin"},
    {"walltime": "08:00:00", "job_name": "npt"},
    {"walltime": "08:00:00", "job_name": "melt"},
    {"walltime": "12:00:00", "job_name": "quench"},
    {"walltime": "48:00:00", "job_name": "nvt"},
]

# default general lammps run resources
GENERAL_QADAPTER = {"walltime": "48:00:00", "job_name": "lmp"}

# default md analysis
ANALYSIS_TYPES = ["diffusion", "rdf"]

# Interface default lammps recipe
INTERFACE_LAMMPS_RECIPE = [
    ["emin", ["template_filename", "min_interface"]],
    ["relax", ["template_filename", "nve_interface"]],
    ["equilibration", ["template_filename", "nvt_interface"]],
    ["production", ["template_filename", "nvt_interface"]],
]

# Interface default lammps recipe settings
INTERFACE_RECIPE_SETTINGS = [
    {
        "data_filename": "../data.mixture",
        "restart_finalname": "restart.emin.restart",
    },
    {
        "restart_filename": "../emin/restart.emin.restart",
        "restart_finalname": "restart.relax.restart"},
    {
        "restart_filename": "../relax/restart.relax.restart",
        "velocity_logic": "",
        "reset_timestep_logic": "",
        "restart_final_filename": "restart.equilibration.restart",
    },
    {
        "restart_filename": "../equilibration/restart.equililibration.restart",
        "reset_timestep_logic": "",
        "restart_final_filename": "restart.production.restart",
        "data_final_filename": "data.production.data",
    },
]

# Interface default lammps run resources
INTERFACE_QADAPTER_RUN_LAMMPS_SPEC = [
    {"walltime": "04:00:00", "job_name": "emin"},
    {"walltime": "24:00:00", "job_name": "rel"},
    {"walltime": "48:00:00", "job_name": "equil"},
    {"walltime": "48:00:00", "job_name": "prod"},
]