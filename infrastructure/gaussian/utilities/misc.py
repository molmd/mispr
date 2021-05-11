# coding: utf-8


# Contains miscellaneous functions useful in many of the infrastructure levels.

import logging

__author__ = 'Rasha Atwi'
__maintainer__ = 'Rasha Atwi'
__email__ = 'rasha.atwi@stonybrook.edu'
__status__ = 'Development'
__date__ = 'Jan 2021'
__version__ = 0.2

logger = logging.getLogger(__name__)


def pass_gout_dict(fw_spec, key):
    gout_dict = fw_spec.get('gaussian_output', {}).get(key)
    proceed_keys = fw_spec.get('proceed', {})
    for k, v in proceed_keys.items():
        if gout_dict['output'].get(k,
                                   gout_dict['output']['output'].get(k)) != v:
            raise ValueError(
                f'The condition for {k} is not met, Terminating'
            )
    return gout_dict


def recursive_signature_remove(d):
    if isinstance(d, dict):
        return {i: recursive_signature_remove(j)
                for i, j in d.items() if not i.startswith('@')}
    else:
        return d


