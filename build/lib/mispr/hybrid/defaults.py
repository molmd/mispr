OPT_GAUSSIAN_INPUTS = {
    "functional": "wB97X",
    "basis_set": "Def2TZVP",
    "route_parameters": {
        "Opt": "(calcfc, tight)",
        "SCF": "Tight",
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
    },
    "link0_parameters": {
        "%chk": "checkpoint.chk",
        "%mem": "45GB",
        "%NProcShared": "28",
    },
}

FREQ_GAUSSIAN_INPUTS = {
    "functional": "wB97X",
    "basis_set": "Def2TZVP",
    "route_parameters": {
        "Freq": None,
        "iop(7/33=1)": None,
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
    },
    "link0_parameters": {
        "%chk": "checkpoint.chk",
        "%mem": "45GB",
        "%NProcShared": "28",
    },
}

NMR_GAUSSIAN_INPUTS = {
    "functional": "wB97X",
    "basis_set": "Def2TZVP",
    "route_parameters": {
        "NMR": "GIAO",
        "iop(7/33=1)": None,
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
    },
    "link0_parameters": {
        "%chk": "checkpoint.chk",
        "%mem": "45GB",
        "%NProcShared": "28",
    },
}
