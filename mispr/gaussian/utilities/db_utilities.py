# coding: utf-8


# Contains db utility functions.

import os
import logging

from fireworks.fw_config import CONFIG_FILE_DIR

from mispr.gaussian.database import GaussianCalcDb

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)


def get_db(input_db=None):
    if not input_db:
        input_db = f"{CONFIG_FILE_DIR}/db.json"
        if not os.path.isfile(input_db):
            raise FileNotFoundError("Please provide the database configurations")
    if isinstance(input_db, dict):
        db = GaussianCalcDb(**input_db)
    else:
        db = GaussianCalcDb.from_db_file(input_db)

    return db


def find_calc_in_db(query_criteria, db):
    from mispr.gaussian.utilities.gout import process_run

    g_out = process_run(operation_type="get_from_run_query", run=query_criteria, db=db)
    pass
