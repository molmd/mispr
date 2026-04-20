import os
import logging

from __future__ import annotations
from typing import Any, Dict, List, Optional

from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from mispr.vasp.defaults import TEMPLATE_TYPES
from mispr.vasp.utilities.structure import process_structure
from mispr.vasp.utilities.poscar import write_poscar
from mispr.vasp.utilities.potcar import generate_potcar

from pymatgen.io.vasp.inputs import Incar, Kpoints, Structure

__author__ = "Sourav Maiti"
__maintainer__ = "Sourav Maiti"
__email__ = "sourav.maiti@stonybrook.edu"
__status__ = "Development"
__date__ = "Mar 2026"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

TEMPLATE_DIR = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "templates")
)
DEFAULT_KEY = "vaspin_key"

@explicit_serialize
class WritePOSCARandPOTCAR(FiretaskBase):
    
    _fw_name = "write POSCAR and POTCAR"
    required_params: List[str] = ["structure"]
    optional_params: List[str] = [
        "working_dir",
        "structure_operation_type",
        "mp_api_key",
        "mp_id",
        "poscar_operation_type",
        "multiplicity",
        "miller_index",
        "min_slab_size",
        "min_vacuum_size",
        "center_slab",
        "primitive",
        "max_normal_search",
        "potcar_dir",
        "potcar_functional",
    ]

    def run_task(self, fw_spec):

        working_dir = fw_spec.get("working_dir",
                                  self.get("working_dir", os.getcwd()))
        os.makedirs(working_dir, exist_ok=True)
        os.chdir(working_dir)

        struct_src = self.get("structure")
        str_op_type = self.get("structure_operation_type", "from_file")
        structure_kwargs = {
            "mp_api_key": self.get("mp_api_key") or os.environ.get("MP_API_KEY"),
            "mp_id": self.get("mp_id"),
        }
        struct = process_structure(str_op_type, struct_src, **structure_kwargs)
        poscar_op_type = self.get("poscar_operation_type", "unitcell")
        poscar_kwargs = {
            "multiplicity": self.get("multiplicity"),
            "miller_index": self.get("miller_index"),
            "min_slab_size": self.get("min_slab_size"),
            "min_vacuum_size": self.get("min_vacuum_size"),
            "center_slab": self.get("center_slab"),
            "primitive": self.get("primitive"),
            "max_normal_search": self.get("max_normal_search"),
        }
        write_poscar(struct, poscar_op_type, **poscar_kwargs)
        potcar_functional = self.get("potcar_functional", "PBE")
        potcar_dir = self.get("potcar_dir", os.environ.get("PMG_VASP_PSP_DIR"))
        generate_potcar("POSCAR", potcar_dir, potcar_functional)

        output_dir = self.get("output_dir") or os.getcwd()

        logger.info(
            "Wrote POSCAR and POTCAR (preset=%s) to %s", self["preset"], output_dir
        )

        return FWAction(
            update_spec={
                "structure": struct.as_dict(),
                "calc_dir": os.path.abspath(output_dir),
            }
        )
    
@explicit_serialize
class WriteINCAR(FiretaskBase):
    _fw_name = "Write INCAR"
    required_params = []
    optional_params = [
        "working_dir",
        "Incar_settings",
        "template_filename",
        "template_dir",
        "template_str",
        "vaspin_key",
    ]

    def run_task(self, fw_spec):
        
        working_dir = fw_spec.get("working_dir", self.get("working_dir",
                                                          os.getcwd()))
        template_dir = None
        os.makedirs(working_dir, exist_ok=True)
        os.chdir(working_dir)
        
        if isinstance(self.get("template_str"), str):
            template_string = self.get("template_str")
            template_filename = None
        elif isinstance(self.get("template_filename"), str):
            template_filename = self.get("template_filename")
            if template_filename in TEMPLATE_TYPES:
                template_dir = TEMPLATE_DIR
            elif isinstance(self.get("template_dir"), str):
                template_dir = self.get("template_dir")
            else:
                raise KeyError(
                    "Directory containing custom template file was not "
                    "specified; add as optional parameter"
                )
        else:
            raise KeyError(
                "Either template was not provided as a valid string or the "
                "path to a template text file was not provided."
            )

        if template_filename:
            template_path = os.path.join(template_dir, template_filename)
            with open(template_path) as file:
                template_string = file.read()   

        incar = Incar(
            template_string,
            )
        
        settings = self.get("Incar_settings", {})
        updates = settings.split(",")

        for item in updates:
            key, value = item.split(":")
            key = key.strip()
            value = value.strip()
            
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    if value.lower() in ["true", "false"]:
                        value = value.lower() == "true"
            incar[key] = value

        incar.write_file("INCAR")

        return FWAction(
            update_spec={
                "calc_dir": os.path.abspath(working_dir),
            }
        )
    
@explicit_serialize
class WriteKPOINTS(FiretaskBase):
    _fw_name = "Write KPOINTS"
    required_params = []
    optional_params = [
        "working_dir",
        "template_filename",
        "template_dir",
        "template_str",
        "kpoints_settings",
        "vaspin_key",
    ]

    def run_task(self, fw_spec):
        
        working_dir = fw_spec.get("working_dir", self.get("working_dir",
                                                          os.getcwd()))
        template_dir = None
        os.makedirs(working_dir, exist_ok=True)
        os.chdir(working_dir)
        
        if isinstance(self.get("template_str"), str):
            template_string = self.get("template_str")
            template_filename = None
        elif isinstance(self.get("template_filename"), str):
            template_filename = self.get("template_filename")
            if template_filename in TEMPLATE_TYPES:
                template_dir = TEMPLATE_DIR
            elif isinstance(self.get("template_dir"), str):
                template_dir = self.get("template_dir")
            else:
                raise KeyError(
                    "Directory containing custom template file was not "
                    "specified; add as optional parameter"
                )
        else:
            raise KeyError(
                "Either template was not provided as a valid string or the "
                "path to a template text file was not provided."
            )

        if template_filename:
            template_path = os.path.join(template_dir, template_filename)
            with open(template_path) as file:
                template_string = file.read()
                kpoints = Kpoints(template_string)
        else:
            kpoints_settings = self.get("kpoints_settings", "automatic, 0, gamma")
            kpoints_settings = kpoints_settings.lower().split()
            if kpoints_settings[0] == "gamma":
                grid = [int(x) for x in kpoints_settings[1:4]]
                kpoints = Kpoints.gamma_automatic(grid)

            elif kpoints_settings[0] == "monkhorst":
                grid = [int(x) for x in kpoints_settings[1:4]]
                kpoints = Kpoints.monkhorst_automatic(grid)
            
            else:
                raise ValueError(
            "Unsupported KPOINTS format. Use: gamma, monkhorst, or provide a template file/string."
        )
        
        kpoints.write_file("KPOINTS")

        return FWAction(
            update_spec={
                "calc_dir": os.path.abspath(working_dir),
            }
        )

