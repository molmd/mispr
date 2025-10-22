"""Define utility functions for handling files and paths."""

import os
import logging

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.4"

logger = logging.getLogger(__name__)


def bibtex_parser(bib_file, working_dir):
    """
    Parse a bibtex file and returns a dictionary of the entries.

    Args:
        bib_file (str): Relative or absolute path to the bibtex file.
        working_dir (str): Name of the working directory where the bibtex file is
            located if bib_file path is relative; else None.

    Returns:
        dict: Dictionary of the entries in the bibtex file.
    """
    try:
        import bibtexparser
    except ModuleNotFoundError:
        raise ImportError(
            "Defining standard electrode potential "
            "references requires bibtexparser to be "
            "installed."
        )
    bib_file = recursive_relative_to_absolute_path(bib_file, working_dir)
    print(bib_file)
    with open(bib_file) as bibfile:
        bp = bibtexparser.load(bibfile)
        entry = bp.entries[0]
        return entry


def recursive_relative_to_absolute_path(operand, working_dir):
    """
    Convert recursively relative file paths to absolute file paths.

    Args:
        operand (str, list, dict): File, list of files, or a dictionary where the values
            are the files; the file(s) path can be relative or absolute.
        working_dir (str): Name of the working directory where the file(s) is/are
            located if operand path is relative; else None.

    Returns:
        str or list or dict: File, list of files, or dict where the values are the
            absolute file paths.
    """
    if isinstance(operand, str):
        if os.path.isabs(operand):
            return operand
        elif os.path.isfile(operand):
            return os.path.join(os.getcwd(), operand)
        else:
            full_path = os.path.join(working_dir, operand)
            if os.path.isfile(full_path):
                return full_path
            else:
                return operand
    elif isinstance(operand, dict):
        return {
            i: recursive_relative_to_absolute_path(j, working_dir)
            for i, j in operand.items()
        }
    elif isinstance(operand, list):
        return [recursive_relative_to_absolute_path(i, working_dir) for i in operand]
    else:
        return operand
