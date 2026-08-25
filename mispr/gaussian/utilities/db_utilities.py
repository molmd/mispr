"""Define db utility functions."""

import os
import logging

from fireworks.fw_config import CONFIG_FILE_DIR

from mispr.gaussian.database import GaussianCalcDb

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.4"

logger = logging.getLogger(__name__)


def get_db(input_db=None):
    """
    Helper function to create a GaussianCalcDb instance from a file or a dict.

    Args:
        input_db (str or dict, optional): Path to db file or a dict containing db info.

    Returns:
        GaussianCalcDb.
    """
    if not input_db:
        input_db = f"{CONFIG_FILE_DIR}/db.json"
        if not os.path.isfile(input_db):
            raise FileNotFoundError("Please provide the database configurations")
    if isinstance(input_db, dict):
        db = GaussianCalcDb(**input_db)
    else:
        db = GaussianCalcDb.from_db_file(input_db)

    return db
