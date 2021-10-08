import os

# default rdf settings
RDF_SETTINGS = {
    "rdf_type": "atomic",
    "r_cut": 20,
    "bin_size": [0.05],
    "filename": os.path.abspath(os.path.join("../../../lammps", "nvt", "log.lammps")),
    "path_or_buff": "rdf.csv",
    "save_mode": True,
}

# default cn settings
CN_SETTINGS = {
    "cn_type": "atomic",
    "bin_size": [0.05],
    "filename": os.path.abspath(os.path.join("../../../lammps", "nvt", "log.lammps")),
    "path_or_buff": "cn.csv",
    "save_mode": True,
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
DIFF_SETTINGS = {"initial_time": None, "final_time": None, "dimension": 3, "save": True,
                 "plot": True, "diff_dist": False}

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
TEMPLATE_TYPES = ["emin_general", "emin_gaff", "npt", "nvt"]

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
    "pair_style": "lj/cut/coul/long",
    "pair_style_args": 10.0,
    "kspace_style": "pppm",
    "kspace_style_args": "1.0e-4",
    "data_file_name": "complex.data",
    "restart_final_name": "restart.emin.restart",
    "dump_file_name": "dump.emin.dump",
    "dump_modify_args": "",
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
    "pair_style_args": 10.0,
    "kspace_style": "pppm",
    "kspace_style_args": "1.0e-4",
    "restart_filename": "restart.emin.restart",
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
    "dump_file_name": "dump.npt.*.dump",
    "dump_modify_logic": "#",
    "dump_modify_args": "",
    "restart_period": 1000000,
    "restart_intermediate_file_name": "restart.npt.*.restart",
    "run": 2000000,
    "restart_final_file_name": "restart.npt.restart",
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
    "pair_style_args": 10.0,
    "kspace_style": "pppm",
    "kspace_style_args": "1.0e-4",
    "restart_file_name": "restart.npt.restart",
    "group_definitions": "",
    "shake_logic": "#",
    "shake_group": "",
    "shake_topologies": "",
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
    "dump_file_name": "dump.nvt.*.dump",
    "dump_modify_logic": "#",
    "dump_modify_args": "",
    "restart_period": 1000000,
    "restart_intermediate_file_name": "restart.nvt.*.restart",
    "run": 5000000,
    "restart_final_file_name": "restart.nvt.restart",
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
    {"data_file_name": "../complex.data", "restart_final_name": "restart.emin.restart"},
    {
        "restart_filename": "../emin/restart.emin.restart",
        "restart_final_file_name": "restart.npt.restart",
    },
    {
        "restart_filename": "../npt/restart.npt.restart",
        "temperature_initial": 500.0,
        "temperature_final": 500.0,
        "run": 1000000,
        "restart_final_file_name": "restart.melt_500K.restart",
    },
    {
        "restart_file_name": "../melt/restart.melt_500K.restart",
        "temperature_initial": 500.0,
        "temperature_final": 298.15,
        "run": 2000000,
        "restart_final_file_name": "restart.quench_298-15K.restart",
    },
    {
        "restart_file_name": "../quench/restart.quench_298-15K.restart",
        "temperature_initial": 298.15,
        "temperature_final": 298.15,
        "run": 5000000,
        "restart_final_file_name": "restart.nvt_5-ns.restart",
    },
]

# default lammps run resources
QADAPTER_RUN_LAMMPS_SPEC = [
    {"walltime": "00:10:00", "job_name": "emin"},
    {"walltime": "08:00:00", "job_name": "npt"},
    {"walltime": "03:00:00", "job_name": "melt"},
    {"walltime": "06:00:00", "job_name": "quench"},
    {"walltime": "15:00:00", "job_name": "nvt"},
]

# default general lammps run resources
GENERAL_QADAPTER = {"walltime": "48:00:00", "job_name": "lmp"}

# default md analysis
ANALYSIS_TYPES = ["diffusion", "rdf"]
