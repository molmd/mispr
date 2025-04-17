# gaussian jobs used to identify the type of calculation and insert it to the db
JOB_TYPES = {
    "sp",
    "opt",
    "freq",
    "irc",
    "ircmax",
    "scan",
    "polar",
    "admp",
    "bomd",
    "eet",
    "force",
    "stable",
    "volume",
    "density",
    "guess",
    "pop",
    "scrf",
    "cphf",
    "prop",
    "nmr",
    "cis",
    "zindo",
    "td",
    "eom",
    "sac-ci",
}

# gaussian SCRF models used to identify calculation model and insert to the db
SCRF_MODELS = {"pcm", "iefpcm", "cpcm", "dipole", "ipcm", "isodensity", "scipcm", "smd"}

# default gaussian inputs used in the workflows if not specified by user
STANDARD_OPT_GUASSIAN_INPUT = {
    "functional": "B3LYP",
    "basis_set": "6-31G(d)",
    "route_parameters": {"Opt": None},
    "link0_parameters": {
        "%chk": "checkpoint.chk",
        "%mem": "45GB",
        "%NProcShared": "24",
    },
}

# maximum number of errors to correct
CUSTODIAN_MAX_ERRORS = 5
