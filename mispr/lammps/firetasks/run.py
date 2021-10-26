# coding: utf-8


# Defines firetasks for running LAMMPS simulations and AmberTools.

import os
import re
import shutil
import logging
import subprocess

from configparser import ConfigParser

from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from pymatgen.io.lammps.inputs import LammpsInput

from mispr.gaussian.utilities.misc import recursive_compare_dicts

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

CONFIG_PATH = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "config", "config.ini"
    )
)


@explicit_serialize
class RunLammpsDirect(FiretaskBase):
    _fw_name = "Run Lammps"
    required_params = []
    optional_params = ["working_dir", "control_filename", "lammps_cmd"]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        control_filename = self.get("control_filename", "complex.lammpsin")
        control_file_path = os.path.join(working_dir, control_filename)

        ntasks_node = fw_spec.get("_queueadapter", {"ntasks_per_node": 1}).get(
            "ntasks_per_node", 1
        )
        nodes = fw_spec.get("_queueadapter", {"nodes": 1}).get("nodes", 1)
        net_ntasks = ntasks_node * nodes
        command = self.get("lammps_cmd")

        if not command:
            config = ConfigParser()
            config.read(CONFIG_FILE_DIR + "/config.ini")
            command = config["LammpsRunCalc"]["lcmd"]
        command = command.replace("$control_path$", control_file_path)
        command = command.replace("$SLURM_NTASKS", str(net_ntasks))

        logger.info("Running command: {}".format(command))
        return_code = subprocess.call(command, shell=True)
        logger.info("Finished running with return code: {}".format(return_code))


@explicit_serialize
class RunLammpsFake(FiretaskBase):
    required_params = ["ref_dir"]
    optional_params = ["working_dir", "control_filename"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    def _verify_inputs(self):
        ref_dir = self["ref_dir"]
        working_dir = self.get("working_dir", os.getcwd())
        control_file = self.get("control_filename", "complex.lammpsin")

        user_control = LammpsInput.from_file(f"{working_dir}/{control_file}")
        ref_control = LammpsInput.from_file(f"{ref_dir}/{control_file}")

        ref_dict = ref_control.as_dict()
        user_dict = user_control.as_dict()

        # remove spaces and keys starting with a #
        user_dict = {
            k: user_dict[k]
            for k in user_dict
            if not re.match(r"#(\s+)?", k)
            if re.match(r"\S+?", k)
        }
        ref_dict = {
            k: ref_dict[k]
            for k in ref_dict
            if not re.match(r"#(\s+)?", k)
            if re.match(r"\S+?", k)
        }
        diff = recursive_compare_dicts(ref_dict, user_dict, "ref_dict", "user_dict")

        if diff:
            raise ValueError(
                f"Control settings are inconsistent with reference control!\n{diff}!"
            )
        logger.info("RunLammpsFake: verified control successfully")

    def _clear_inputs(self):
        working_dir = self.get("working_dir", os.getcwd())
        control_file = self.get("control_filename", "complex.lammpsin")
        control_file = f"{working_dir}/{control_file}"
        if os.path.exists(control_file):
            os.remove(control_file)

    def _generate_outputs(self):
        ref_dir = self["ref_dir"]
        working_dir = self.get("working_dir", os.getcwd())
        for file in os.listdir(ref_dir):
            full_path = f"{ref_dir}/{file}"
            if os.path.isfile(full_path):
                shutil.copy(full_path, working_dir)
        logger.info("RunLammpsFake: ran fake Lammps, generated outputs")


# TODO: LAMMPS Custodian Firetask


@explicit_serialize
class RunAntechamber(FiretaskBase):
    _fw_name = "Run Antechamber"
    required_params = []
    optional_params = [
        "working_dir",
        "input_filename_a",
        "input_file_type",
        "output_filename_a",
        "output_file_type",
        "charge_method",
        "antechamber_cmd",
    ]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        input_filename = self.get("input_filename_a", "mol.esp")
        input_file_type = self.get("input_filetype", "gesp")
        output_filename = self.get("output_filename_a", "mol.mol2")
        output_file_type = self.get("output_filetype", "mol2")
        charge_method = self.get("charge_method", "resp")

        input_file_path = os.path.join(working_dir, input_filename)
        output_file_path = os.path.join(working_dir, output_filename)

        command = self.get("antechamber_cmd")

        if not command:
            config = ConfigParser()
            config.read(CONFIG_FILE_DIR + "/config.ini")
            command = config["AmbertoolsRunCalc"]["acmd"]
        command = (
            command.replace("$input_file$", input_file_path)
            .replace("$input_type$", input_file_type)
            .replace("$output_file$", output_file_path)
            .replace("$output_type$", output_file_type)
            .replace("$charge_method$", charge_method)
        )

        logger.info("Running command: {}".format(command))
        return_code = subprocess.call(command, shell=True)
        logger.info("Finished running with return code: {}".format(return_code))


# TODO: Antechamber Custodian Firetask


@explicit_serialize
class RunParmchk(FiretaskBase):
    _fw_name = "Run Parmchk"
    required_params = []
    optional_params = [
        "working_dir",
        "input_filename_p",
        "output_filename_p",
        "parmchk_cmd",
    ]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        input_filename = self.get("input_filename_p", "mol.mol2")
        output_filename = self.get("output_filename_p", "mol.frcmod")

        input_file_path = os.path.join(working_dir, input_filename)
        output_file_path = os.path.join(working_dir, output_filename)

        command = self.get("parmchk_cmd")

        if not command:
            config = ConfigParser()
            config.read(CONFIG_FILE_DIR + "/config.ini")
            command = config["AmbertoolsRunCalc"]["pcmd"]
        command = command.replace("$input_file$", input_file_path).replace(
            "$output_file$", output_file_path
        )

        logger.info("Running command: {}".format(command))
        return_code = subprocess.call(command, shell=True)
        logger.info("Finished running with return code: {}".format(return_code))


# TODO: Edit RunParmchk for general input file type
# TODO: Parmchk Custodian Firetask


@explicit_serialize
class RunTleap(FiretaskBase):
    _fw_name = "Run Tleap"
    required_params = []
    optional_params = ["working_dir", "script_filename", "tleap_cmd"]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        script_filename = self.get("script_filename", "tleap.in")
        script_file_path = os.path.join(working_dir, script_filename)

        command = self.get("tleap_cmd")

        if not command:
            config = ConfigParser()
            config.read(CONFIG_FILE_DIR + "/config.ini")
            command = config["AmbertoolsRunCalc"]["tcmd"]
        command = command.replace("$input_file$", script_file_path)

        logger.info("Running command: {}".format(command))
        return_code = subprocess.call(command, shell=True)
        logger.info("Finished running with return code: {}".format(return_code))


# TODO: tleap Custodian Firetask
