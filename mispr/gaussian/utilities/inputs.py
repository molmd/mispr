# coding: utf-8


# Contains functions for handling gaussian inputs.

import logging

from copy import deepcopy

from mispr.gaussian.defaults import STANDARD_OPT_GUASSIAN_INPUT

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)


def _add_solvent_inputs(
    gaussian_inputs, solvent_gaussian_inputs, solvent_properties=None
):
    if "generic" in solvent_gaussian_inputs.lower() and not solvent_properties:
        raise Exception(
            "A generic solvent is provided as an input without "
            "specifying its parameters."
        )
    for key, value in gaussian_inputs.items():
        if "route_parameters" not in value:
            gaussian_inputs[key]["route_parameters"] = {}
        gaussian_inputs[key]["route_parameters"]["SCRF"] = solvent_gaussian_inputs
        if solvent_properties:
            if "input_parameters" not in gaussian_inputs[key]:
                gaussian_inputs[key]["input_parameters"] = {}
            gaussian_inputs[key]["input_parameters"].update(solvent_properties)
    return gaussian_inputs


def _check_solvent_inputs(gaussian_inputs):
    # gaussian_inputs is a list of dicts
    route_params = {}
    for key, value in gaussian_inputs.items():
        if value:
            route_params.update(value.get("route_parameters", {}))
    assert "scrf" not in [i.lower() for i in route_params], (
        "solvent inputs should be provided as separate inputs via "
        "solvent_gaussian_inputs and solvent_properties"
    )


def _get_gaussian_inputs(gaussian_inputs, supported_jobs=None):
    # gaussian_inputs is a dict of dicts: {'opt': {}, 'freq': {}, 'nmr': {},
    # 'esp': {}, 'sp': {}}; this function is meant to be used in workflows in
    # which multiple Gaussian jobs are performed and the jobs share similar
    # Gaussian keywords; used to handle situations in which the user is not
    # explicitly defining every single job input dictionary
    supported_jobs = supported_jobs or {}
    supported_jobs = {
        **{
            "freq": {"Freq": None},
            "nmr": {"NMR": "GIAO"},
            "esp": {"pop": "MK", "iop(6/50=1)": None},
            "sp": {"SP": None},
        },
        **{k.lower(): v for k, v in supported_jobs.items()},
    }

    gaussian_inputs = {
        k.lower(): v if v is not None else {} for k, v in gaussian_inputs.items()
    }

    if "opt" not in gaussian_inputs:
        gaussian_inputs["opt"] = {}
    gaussian_inputs["opt"] = {**STANDARD_OPT_GUASSIAN_INPUT, **gaussian_inputs["opt"]}
    if "opt" not in [i.lower() for i in gaussian_inputs["opt"]["route_parameters"]]:
        gaussian_inputs["opt"]["route_parameters"].update({"Opt": None})

    for job in gaussian_inputs:
        if job in supported_jobs and job != "opt":
            gaussian_inputs[job] = _update_gaussian_inputs(
                deepcopy(gaussian_inputs["opt"]),
                gaussian_inputs[job],
                supported_jobs[job],
            )
        elif job == "opt":
            pass
        else:
            logger.error(
                "keyword for {} is not known. Please add keyword to "
                "the supported_jobs dict.".format(job)
            )
    return gaussian_inputs


def _update_gaussian_inputs(opt_gaussian_inputs, other_gaussian_inputs, main_keyword):
    gaussian_inputs = {**opt_gaussian_inputs, **other_gaussian_inputs}
    if list(main_keyword.keys())[0].lower() not in [
        i.lower() for i in gaussian_inputs["route_parameters"]
    ]:
        gaussian_inputs["route_parameters"].update(main_keyword)
    for i in gaussian_inputs["route_parameters"]:
        if i.lower() == "opt":
            del gaussian_inputs["route_parameters"][i]
            break
    return gaussian_inputs


def handle_gaussian_inputs(
    gaussian_inputs, solvent_gaussian_inputs=None, solvent_properties=None
):
    gaussian_inputs = _get_gaussian_inputs(gaussian_inputs)
    _check_solvent_inputs(gaussian_inputs)
    if solvent_gaussian_inputs:
        gaussian_inputs = _add_solvent_inputs(
            gaussian_inputs, solvent_gaussian_inputs, solvent_properties
        )
    return gaussian_inputs
