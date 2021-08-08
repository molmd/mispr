import os

# default rdf settings
RDF_SETTINGS = {
    "r_cut": 20,
    "bin_size": [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
    "filename": os.path.abspath(os.path.join("../../../lammps", "nvt", "log.lammps")),
    "path_or_buff": "rdf.csv",
    "save_mode": True,
}

# default msd settings
MSD_SETTINGS = {
    "dt": 1,
    "tao_coeff": 4,
    "msd_type": "com",
    "return_all": False,
    "com_drift": False,
    "avg_interval": False,
    "save_msd": True,
}

# default diffusion settings
DIFF_SETTINGS = {"initial_time": None, "final_time": None, "dimension": 3}
