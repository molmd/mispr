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
__version__ = "0.0.3"

logger = logging.getLogger(__name__)


def _add_solvent_inputs(
    gaussian_inputs, solvent_gaussian_inputs, solvent_properties=None
):
    """
    Adds solvent inputs to the gaussian_inputs dict. The input
    parameters relevant to implicit solvation are required to be
    provided separate from the gaussian_inputs dict to the Gaussian
    workflows. This is done in order to easily save the solvent
    properties to the database or to the json files generated from the
    workflows.

    Args:
        gaussian_inputs (dict): dictionary of dictionaries of Gaussian
            inputs for one or more Gaussian job type, e.g.
            {"opt": {
                "functional": "B3LYP",
                "basis_set": "6-31G",
                "route_parameters": {
                    "Opt": "(calcfc, tight)",
                    "test": None},
                    }
            }
        solvent_gaussian_inputs (str): string of Gaussian inputs for
            the solvent, e.g. "(Solvent=Generic, Read)"
        solvent_properties (dict): dictionary of solvent properties,
            e.g. {"Eps": 4.33, "EpsInf": 1.69}

    Returns:
        dict: dictionary of dictionaries of Gaussian inputs with
            solvent inputs added
    """
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
    """
    Ensures that implicit solvent parameters are not included in the
    Gaussian input dictionaries, since those should be provided
    separately.

    Args:
        gaussian_inputs (dict): dictionary of dictionaries of Gaussian
            inputs for one or more Gaussian job type, e.g.
            {"opt": {
                "functional": "B3LYP",
                "basis_set": "6-31G",
                "route_parameters": {
                    "Opt": "(calcfc, tight)",
                    "test": None},
                    }
            }
    """
    route_params = {}
    for key, value in gaussian_inputs.items():
        if value:
            route_params.update(value.get("route_parameters", {}))
    assert "scrf" not in [i.lower() for i in route_params], (
        "solvent inputs should be provided as separate inputs via "
        "solvent_gaussian_inputs and solvent_properties"
    )


def _get_gaussian_inputs(gaussian_inputs, supported_jobs=None):
    """
    This function is meant to be used in workflows in which multiple
    Gaussian jobs are performed and the jobs share similar Gaussian
    keywords; used to handle situations in which the user is
    not explicitly defining every single job input dictionary.

    Args:
        gaussian_inputs (dict): dictionary of dictionaries of Gaussian
            input parameters for different job types in a given workflow,
            e.g. {"opt": {}, "freq": {}, "nmr": {}, "esp": {}, "sp": {}}
        supported_jobs (dict): dictionary of dictionaries of supported
            job types and their main inputs,
            e.g. {"opt": {"Opt": None}, "nmr": {"NMR": "GIAO"}}

    Returns:
        dict: dictionary of dictionaries of Gaussian input parameters
    """
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
    """
    Uses the fully defined optimization input parameters to fill in the
    input parameters for other job types in a given workflow. Done to
    avoid having to define every single job input dictionary.

    Args:
        opt_gaussian_inputs (dict): dictionary of Gaussian input
            parameters for the optimization step of a given workflow
        other_gaussian_inputs (dict): dictionary of Gaussian inputs
            other than the route_parameters (e.g. input_parameters) for
            a job other than optimization in a given workflow
        main_keyword (dict): simple dictionary containing the main
            Gaussian route_parameters for the job type,
            e.g. {"Freq": None} for a frequency analysis

    Returns:

    """
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
    """
    Wrapper function to cleanup/modify the Gaussian input parameters
    for one or more job in a workflow. Checks for implicit solvent
    parameters and adds missing keywords for a given job.

    Args:
        gaussian_inputs (dict): dictionary of dictionaries of Gaussian
            inputs, e.g. {"opt": {opt_gaussian_inputs},
                          "freq": {freq_gaussian_inputs},
                          etc.}
        solvent_gaussian_inputs (str): string of Gaussian inputs for
            the solvent, e.g. "(Solvent=Generic, Read)"
        solvent_properties (dict): dictionary of solvent properties,
            e.g. {"Eps": 4.33, "EpsInf": 1.69}

    Returns:
        dict: dictionary of dictionaries of reformatted Gaussian inputs
    """
    gaussian_inputs = _get_gaussian_inputs(gaussian_inputs)
    _check_solvent_inputs(gaussian_inputs)
    if solvent_gaussian_inputs:
        gaussian_inputs = _add_solvent_inputs(
            gaussian_inputs, solvent_gaussian_inputs, solvent_properties
        )
    return gaussian_inputs
