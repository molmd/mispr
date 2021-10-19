# coding: utf-8


# Contains miscellaneous functions useful in many of the mispr levels.

import logging

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)


def pass_gout_dict(fw_spec, key):
    gout_dict = fw_spec.get("gaussian_output", {}).get(key)
    proceed_keys = fw_spec.get("proceed", {})
    for k, v in proceed_keys.items():
        if gout_dict["output"].get(k, gout_dict["output"]["output"].get(k)) != v:
            raise ValueError(f"The condition for {k} is not met, Terminating")
    return gout_dict


def recursive_signature_remove(d):
    if isinstance(d, dict):
        return {
            i: recursive_signature_remove(j)
            for i, j in d.items()
            if not i.startswith("@")
        }
    else:
        return d


def recursive_compare_dicts(dict1, dict2, dict1_name, dict2_name, path=""):
    error = ""
    old_path = path
    for key in dict1.keys():
        path = f"{old_path}[{key}]"
        if key not in dict2.keys():
            error += f"Key {dict1_name}{path} not in {dict2_name}\n"
        else:
            if isinstance(dict1[key], dict) and isinstance(dict2[key], dict):
                error += recursive_compare_dicts(
                    dict1[key], dict2[key], "d1", "d2", path
                )
            else:
                if dict1[key] != dict2[key]:
                    error += (
                        f"Value of {dict1_name}{path} ({dict1[key]}) "
                        f"not same as {dict2_name}{path} ({dict2[key]})\n"
                    )

    for key in dict2.keys():
        path = f"{old_path}[{key}]"
        if key not in dict1.keys():
            error += f"Key {dict2_name}{path} not in {dict1_name}\n"
    return error
