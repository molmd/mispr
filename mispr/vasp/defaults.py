import os

# vasp recipt template files
# TODO: decide if I want specific templates for INCAR and KPOINTS
#TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")
TEMPLATE_TYPES = ["INCAR", "KPOINTS"]

# default emin settings
INCAR_SETTINGS = {
    "SYSTEM": "system",
    "ISTART": 0,
    "NCORE": 1,
    "ISMEAR": -5,
    "SIGMA": 0.05,
    "IBRION": 2,
    "ICHARG": 2,
    "NSW": 99,
    "NELM": 100,
    "PREC": "Accurate",
    "ALGO": "Fast",
    "GGA": "PE",
    "ISIF": 3,
    "ENCUT": 520,
    "EDIFF": "1E-5",
    "ISPIN": 2,
    "LREAL": "A",
    "KPOINTS_BSE": "-1,0,0,0",
    "LWAVE": "F",
    "LCHARG": "F",
    "LORBIT": 11,
}

#deafult kpoint settings
KPOINTS_SETTINGS = {
    "gamma": {"kpoints_type": "Gamma", "kpoints_grid": [1, 1, 1]},
    "monkhorst_pack": {"kpoints_type": "Monkhorst-Pack", "kpoints_grid": [3, 3, 3]},
}

# default vasp recipe
VASP_RECIPE = [
    ["INCAR", ["template_filename", "INCAR"]],
    ["KPOINTS", ["template_filename", "KPOINTS"]],
]

# default vasp recipe settings
RECIPE_SETTINGS = [
    {"POSCAR": "POSCAR", "POTCAR": "POTCAR"},
]

# default lammps run resources
QADAPTER_RUN_LAMMPS_SPEC = [
    {"walltime": "48:00:00", "job_name": "Optimize"},
]

# default general lammps run resources
GENERAL_QADAPTER = {"walltime": "48:00:00", "job_name": "vasp_run"}
