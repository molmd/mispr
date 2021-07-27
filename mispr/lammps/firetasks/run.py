# coding: utf-8

# Defines firetasks for running LAMMPS simulations and AmberTools

import os
import logging
import subprocess

from configparser import ConfigParser

from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@tufts.edu"
__status__ = "Development"
__date__ = "Apr 14, 2020"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

CONFIG_PATH = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "config", "config.ini"
    )
)


@explicit_serialize
class RunLammps(FiretaskBase):
    _fw_name = "Run Lammps"
    required_params = []
    optional_params = ["working_dir", "control_filename", "lammps_cmd"]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        control_filename = self.get("control_filename", "complex.lammpsin")
        control_file_path = os.path.join(working_dir, control_filename)

        command = self.get("lammps_cmd")

        if not command:
            config = ConfigParser()
            config.read(CONFIG_PATH)
            command = config["LammpsRunCalc"]["lcmd"]
        command = command.replace("$control_path$", control_file_path)

        logger.info("Running command: {}".format(command))
        return_code = subprocess.call(command, shell=True)
        logger.info("Finished running with return code: {}".format(return_code))


# TODO: LAMMPS Custodian Firetask

# TODO: Fake LAMMPS Firetask


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
            config.read(CONFIG_PATH)
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

        # if isinstance(fw_spec.get("tleap_settings", {}), dict):
        #     tleap_settings = fw_spec.get("tleap_settings", {}).update({})
        #     return FWAction()


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
            config.read(CONFIG_PATH)
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
            config.read(CONFIG_PATH)
            command = config["AmbertoolsRunCalc"]["tcmd"]
        command = command.replace("$input_file$", script_file_path)

        logger.info("Running command: {}".format(command))
        return_code = subprocess.call(command, shell=True)
        logger.info("Finished running with return code: {}".format(return_code))


# TODO: tleap Custodian Firetask
