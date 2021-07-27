# coding: utf-8


# Contains functions for processing different gaussian output formats.

import os
import json
import logging

from bson import ObjectId

from pymatgen.io.gaussian import GaussianOutput

from mispr.gaussian.utilities.misc import recursive_signature_remove
from mispr.gaussian.utilities.dbdoc import _cleanup_gout
from mispr.gaussian.utilities.db_utilities import get_db

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)


def process_run(operation_type, run, input_file=None, **kwargs):
    working_dir = kwargs.get("working_dir", os.getcwd())

    def get_db_():
        return get_db(kwargs["db"]) if "db" in kwargs else get_db()

    if operation_type == "get_from_gout":
        if not isinstance(run, GaussianOutput):
            raise Exception(
                "run is not a GaussianOutput object; either "
                "provide a GaussianOutput object or use another "
                "operation type with its corresponding inputs"
            )
        gout = run.as_dict()
        gout_dict = _cleanup_gout(gout, working_dir, input_file)

    elif operation_type == "get_from_gout_file":
        if not os.path.isabs(run):
            file_path = os.path.join(working_dir, run)
        else:
            file_path = run
        if not os.path.exists(file_path):
            raise Exception(
                "run is not a valid path; either provide a valid "
                "path or use another operation type with its "
                "corresponding inputs"
            )
        try:
            gout = GaussianOutput(file_path).as_dict()
            gout_dict = _cleanup_gout(gout, working_dir, input_file)
        except IndexError:
            raise ValueError(
                "run is not a Gaussian output file; either "
                "provide a valid Gaussian output file or use "
                "another operation type with its corresponding "
                "inputs"
            )

    elif operation_type == "get_from_run_dict":
        if not isinstance(run, dict) and "output" not in run:
            raise Exception(
                "run is not a GaussianOutput dictionary; either"
                "provide a GaussianOutput dictionary or use another"
                "operation type with its corresponding inputs"
            )
        gout_dict = run

    elif operation_type == "get_from_run_id":
        # run = run_id
        db = get_db_()
        run = db.runs.find_one({"_id": ObjectId(run)})
        if not run:
            raise Exception("Gaussian run is not in the database")
        gout_dict = run

    elif operation_type == "get_from_run_query":
        # run = {'inchi': inchi, 'smiles': smiles, 'type': type,
        #        'functional': func, 'basis': basis, 'phase': phase, ...}
        logger.info(
            "If the query criteria satisfy more than "
            "one document, the last updated one will "
            "be used. To perform a more specific "
            "search, provide the document id using "
            "gout_id"
        )
        db = get_db_()
        run = db.retrieve_run(**run)
        if not run:
            raise Exception("Gaussian run is not in the database")
        run = max(run, key=lambda i: i["last_updated"])
        gout_dict = run

    else:
        raise ValueError(f"operation type {operation_type} is not supported")
    if "_id" in gout_dict:
        gout_dict["_id"] = str(gout_dict["_id"])
    if "last_updated" in gout_dict:
        del gout_dict["last_updated"]
    gout_dict = json.loads(json.dumps(gout_dict))
    gout_dict = recursive_signature_remove(gout_dict)

    return gout_dict
