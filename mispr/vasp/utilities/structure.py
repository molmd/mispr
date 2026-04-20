from __future__ import annotations

import os
import logging

from pymatgen.core import Structure

logger = logging.getLogger(__name__)

def process_structure(operation_type, structure, **kwargs):
    working_dir = kwargs.get("working_dir", os.getcwd())
    
    if operation_type == "from_file":
        if not os.path.isabs(structure):
            file_path = os.path.join(working_dir, structure)
        else:
            file_path = structure
        if not os.path.exists(file_path):
            raise Exception(
                "structure is not a valid path; either provide a valid "
                "path or use another operation type with its "
                "corresponding inputs"
            )
        output_structure = Structure.from_file(file_path)
        
    elif operation_type == "from_mp":
        key = kwargs.get("api_key") or os.environ.get("MP_API_KEY")
        mp_id = structure or kwargs.get("mp_id")

        if mp_id is None:
            raise ValueError(
                "No Materials Project ID provided."
                )
        
        if key is None:
            raise EnvironmentError(
                "No Materials Project API key found.  "
                "Please provide an API key via the 'api_key' argument or set the 'MP_API_KEY' environment variable."
                )
        
        try:
            from mp_api.client import MPRester
            with MPRester(key) as mpr:
                struc = mpr.get_structure_by_material_id(mp_id)
        except ImportError:
            from pymatgen.ext.matproj import MPRester
            with MPRester(key) as mpr:
                struc = mpr.get_structure_by_material_id(mp_id)
        output_structure = Structure.from_dict(struc)        
    
    else:
        raise ValueError(f"Unknown operation_type '{operation_type}'.")
    
    return output_structure

