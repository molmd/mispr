import os
import re
import json
import shutil
import logging
import subprocess

from configparser import ConfigParser

from pymatgen.io.vasp.inputs import Incar, Kpoints

from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from mispr.lammps.utilities.opls import MaestroRunner

__author__ = "Sourav Maiti"
__maintainer__ = "Sourav Maiti"
__email__ = "sourav.maiti@stonybrook.edu"
__status__ = "Development"
__date__ = "Mar 2026"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

CONFIG_PATH = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "config", "config.ini"
    )
)

@explicit_serialize
class RunVASP(FiretaskBase):

    _fw_name = "Run VASP"
    required_params = []
    optional_params = [
        "working_dir",  
        "vasp_cmd", 
        "net_ntasks"
    ]

    def run_task(self, fw_spec):
        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        ntasks_node = fw_spec.get("_queueadapter", {"ntasks_per_node": 1}).get(
            "ntasks_per_node", 1
        )
        nodes = fw_spec.get("_queueadapter", {"nodes": 1}).get("nodes", 1)
        net_ntasks = self.get("net_ntasks", ntasks_node * nodes)
        command = self.get("vasp_cmd")

        if not command:
            config = ConfigParser()
            config.read(CONFIG_FILE_DIR + "/config.ini")
            command = config["VASPRunCalc"]["vcmd"]
        command = command.replace("$SLURM_NTASKS", str(net_ntasks))

        logger.info("Running command: {}".format(command))
        return_code = subprocess.call(command, shell=True)
        logger.info(
            "Finished running with return code: {}".format(return_code)
        )
        return FWAction(
            update_spec={
                "vasp_return_code": return_code,
                "vasp_command": command,
            }
        )
