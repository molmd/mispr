# coding: utf-8


# Contains functions for cleaning up JSON documents.

import os
import logging

from pymatgen.io.gaussian import GaussianInput
from pymatgen.core.structure import Molecule

from mispr.gaussian.defaults import JOB_TYPES, SCRF_MODELS
from mispr.gaussian.utilities.metadata import get_chem_schema

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.4"

logger = logging.getLogger(__name__)


def _job_types(gin):
    """
    Determines the type of the Gaussian job (e.g. opt, freq, etc.)
    from a Gaussian input dictionary

    Args:
        gin (dict): Gaussian input dictionary

    Returns:
        list: list of Gaussian job types
    """
    return sorted(
        list(
            filter(
                lambda x: x
                in {k.lower(): v for k, v in gin["route_parameters"].items()},
                JOB_TYPES,
            )
        )
    )


def _modify_gout(gout):
    """
    Modifies the Gaussian output dictionary by removing unnecessary
    keys and restructuring its format to match the schema of the db

    Args:
        gout (dict): Gaussian output dictionary

    Returns:
        dict: modified Gaussian output dictionary

    """
    gout["input"]["charge"] = gout["charge"]
    gout["input"]["spin_multiplicity"] = gout["spin_multiplicity"]
    del_keys_out = (
        "nsites",
        "unit_cell_formula",
        "reduced_cell_formula",
        "pretty_formula",
        "elements",
        "nelements",
        "charge",
        "spin_multiplicity",
    )
    [gout.pop(k, None) for k in del_keys_out]
    return gout


def _create_gin(gout, working_dir, input_file):
    """
    Creates a Gaussian input dictionary to be used for creating
    Gaussian documents in the db; if an input_file is provided,
    content of the dictionary is mostly taken from it; otherwise,
    uses the gout to create the dictionary

    Args:
        gout (dict): Gaussian output dictionary
        working_dir (str): working directory of the Gaussian input_file
        input_file (str): relative or absolute path to the input_file

    Returns:
        dict: Gaussian input dictionary
    """
    if input_file:
        if not os.path.isabs(input_file):
            input_path = os.path.join(working_dir, input_file)
        else:
            input_path = input_file
        gin = GaussianInput.from_file(input_path).as_dict()
        gin["nbasisfunctions"] = gout["input"]["nbasisfunctions"]
        gin["pcm_parameters"] = gout["input"]["pcm_parameters"]
        return gin
    else:
        gin = gout["input"]
        gin["input_parameters"] = None
        gin["@class"] = "GaussianInput"
        gin["@module"] = "pymatgen.io.gaussian"
        logger.info(
            "input parameters at the end of the Gaussian input "
            "section will not be saved to the database due to "
            "a missing input file"
        )
        return gin


def _cleanup_gout(gout, working_dir, input_file):
    """
    Cleans up the Gaussian dictionary to be saved in the database

    Args:
        gout (dict): Gaussian output dictionary
        working_dir (str): working directory of the Gaussian input_file
        input_file (str): relative or absolute path to the input_file

    Returns:
        dict: cleaned up Gaussian dictionary
    """
    gout = _modify_gout(gout)
    gin = _create_gin(gout, working_dir, input_file)
    del gout["input"]
    gauss_version = gout["output"]["gauss_version"]
    del gout["output"]["gauss_version"]
    job_types = _job_types(gin)
    mol = Molecule.from_dict(gout["output"]["molecule"])
    gout_dict = {
        "input": gin,
        "output": gout,
        "functional": gin["functional"],
        "basis": gin["basis_set"],
        "phase": "solution" if gout["is_pcm"] else "gas",
        "type": ";".join(job_types),
        **get_chem_schema(mol),
        "gauss_version": gauss_version,
    }
    gout_dict = {
        i: j
        for i, j in gout_dict.items()
        if i not in ["sites", "@module", "@class", "charge", "spin_multiplicity"]
    }
    return gout_dict


def add_solvent_to_prop_dict(prop_dict, solvent_gaussian_inputs, solvent_properties):
    """
    Adds solvent properties to a property dictionary (e.g. BDE, BE, etc.)

    Args:
        prop_dict (dict): property dictionary
        solvent_gaussian_inputs (str): Gaussian input parameters
            corresponding to the implicit solvent model used in
            the Gaussian calculations, e.g. "(Solvent=TetraHydroFuran)"
        solvent_properties (dict): additional solvent input parameters
            used in the Gaussian calculations; e.g., {"EPS":12}

    Returns:
        dict: property dictionary with solvent properties added
    """
    if not solvent_gaussian_inputs:
        solvent = "water"
        solvent_model = "pcm"
    else:
        solvent_inputs = [
            i.lower() for i in solvent_gaussian_inputs.strip("()").split(",")
        ]
        solvent = [
            string.split("=")[1] for string in solvent_inputs if "solvent" in string
        ] or ["water"]
        solvent_model = list(
            filter(lambda x: x in {i for i in solvent_inputs}, SCRF_MODELS)
        ) or ["pcm"]
    prop_dict["solvent"] = "".join(solvent)
    prop_dict["solvent_model"] = "".join(solvent_model)
    prop_dict["solvent_properties"] = solvent_properties
    return prop_dict
