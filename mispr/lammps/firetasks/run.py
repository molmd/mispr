# coding: utf-8


# Defines firetasks for running LAMMPS simulations and AmberTools.

import os
import re
import json
import shutil
import logging
import subprocess

from configparser import ConfigParser

from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize

from pymatgen.io.lammps.inputs import LammpsInput

from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.utilities.misc import recursive_compare_dicts
from mispr.gaussian.utilities.metadata import get_mol_formula
from mispr.lammps.utilities.opls import MaestroRunner
from mispr.lammps.utilities.utilities import get_db, process_ff_doc, \
    add_ff_labels_to_dict

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

CONFIG_PATH = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 
        "..", "config", "config.ini"
    )
)

OPLS_DOI = "10.1002/jcc.20292"


@explicit_serialize
class RunLammpsDirect(FiretaskBase):
    """
    Run a LAMMPS simulation from a single control file.
    
    Args:
        working_dir (str, optional): The directory to run the simulation
            in. Defaults to the current working directory. file. 
            Defaults to "complex.lammpsin".
        lammps_cmd (str, optional): The command to run LAMMPS. If not 
            provided, the command will attempt to read from from the 
            config.ini file.
        net_ntasks (int, optional): The number of processors to run the 
            simulation on. If not provided, the number of processors 
            will be read from the fw_spec based on the number of nodes 
            used and the number of tasks per node.
    """
    _fw_name = "Run Lammps"
    required_params = []
    optional_params = [
        "working_dir", 
        "control_filename", 
        "lammps_cmd", 
        "net_ntasks"
    ]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        control_filename = self.get("control_filename", "complex.lammpsin")
        control_file_path = os.path.join(working_dir, control_filename)

        ntasks_node = fw_spec.get("_queueadapter", {"ntasks_per_node": 1}).get(
            "ntasks_per_node", 1
        )
        nodes = fw_spec.get("_queueadapter", {"nodes": 1}).get("nodes", 1)
        net_ntasks = self.get("net_ntasks", ntasks_node * nodes)
        command = self.get("lammps_cmd")

        if not command:
            config = ConfigParser()
            config.read(CONFIG_FILE_DIR + "/config.ini")
            command = config["LammpsRunCalc"]["lcmd"]
        command = command.replace("$control_path$", control_file_path)
        command = command.replace("$SLURM_NTASKS", str(net_ntasks))

        logger.info("Running command: {}".format(command))
        return_code = subprocess.call(command, shell=True)
        logger.info(
            "Finished running with return code: {}".format(return_code)
        )


@explicit_serialize
class RunLammpsFake(FiretaskBase):
    """
    Run a fake LAMMPS simulation.

    Args:
        ref_dir (str): The directory containing the reference LAMMPS 
            files.
        working_dir (str, optional): The directory to run the fake 
            simulation in. Defaults to the current working directory.
        control_filename (str, optional): The name of the LAMMPS control
            file. Defaults to "complex.lammpsin".
    """
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
        diff = recursive_compare_dicts(ref_dict, user_dict, 
                                       "ref_dict", "user_dict")

        if diff:
            raise ValueError(
                f"Control settings are inconsistent with reference \
                    control!\n{diff}!"
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
    """
    Runs the Antechamber program. Intented to be used to generate the 
    partial charges for the molecule in a mol2 file from a Gaussian ESP 
    file. Refer to the AmberTools documentation for more information on 
    this program: https://ambermd.org/AmberTools.php

    Args:
        working_dir (str, optional): The directory to run Antechamber 
            in. Defaults to the current working directory.
        input_filename_a (str, optional): The name of the input file. 
            Defaults to "mol.esp".
        input_file_type (str, optional): The type of the input file. 
            Defaults to "gesp".
        output_filename_a (str, optional): The name of the output file. 
            Defaults to "mol.mol2".
        output_file_type (str, optional): The type of the output file. 
            Defaults to "mol2".
        charge_method (str, optional): The method to calculate the 
            partial charges. Defaults to "resp".
        antechamber_cmd (str, optional): The command to run Antechamber.
            If not provided, the command will be read from the 
            config.ini file.
    """
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
        logger.info(
            "Finished running with return code: {}".format(return_code)
        )


# TODO: Antechamber Custodian Firetask


@explicit_serialize
class RunParmchk(FiretaskBase):
    """
    Runs the Parmchk program. Intended to be used to generate a frcmod 
    file from a mol2 file. Refer to the AmberTools documentation for 
    more information on the Parmchk program: 
    https://ambermd.org/AmberTools.php

    Args:
        working_dir (str, optional): The directory to run Parmchk in. 
            Defaults to the current working directory.
        input_filename_p (str, optional): The name of the input file. 
            Defaults to "mol.mol2".
        output_filename_p (str, optional): The name of the output file. 
            Defaults to "mol.frcmod".
        parmchk_cmd (str, optional): The command to run Parmchk. If not 
            provided, the command will be read from the config.ini file.
    """
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
        logger.info(
            "Finished running with return code: {}".format(return_code)
        )


# TODO: Edit RunParmchk for general input file type
# TODO: Parmchk Custodian Firetask


@explicit_serialize
class RunTleap(FiretaskBase):
    """
    Runs the Tleap program. Intended to be used to generate the GAFF 
    parameters for a single molecular species into a prmtop file from a 
    mol2 file and a frcmod file. Refer to the AmberTools documentation
    for more information on this program: 
    https://ambermd.org/AmberTools.php

    Args:
        working_dir (str, optional): The directory to run Tleap in. 
            Defaults to the current working directory.
        script_filename (str, optional): The name of the script file. 
            Defaults to "tleap.in".
        tleap_cmd (str, optional): The command to run Tleap. If not 
            provided, the command will be read from the config.ini file.
    """
    _fw_name = "Run Tleap"
    required_params = []
    optional_params = [
        "working_dir",
        "script_filename",
        "tleap_cmd",
        "system_force_field_dict",
    ]

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
        logger.info(
            "Finished running with return code: {}".format(return_code)
        )


# TODO: evaluate if the "system_force_field_dict" is necessary

@explicit_serialize
class RunMaestro(FiretaskBase):
    """
    Run the Maestro program using the MaestroRunner class in the 
    mispr.lammps.utilities.opls module to generate OPLS force field 
    parameters for a molecule. After generating the parameters, add 
    these parameters to the spec. Refer to the Maestro documentation for 
    more information on this program: 
    https://www.schrodinger.com/maestro

    Args:
        input_file (str): The path to the input file.
        label (str, optional): The label for the molecule. Defaults to 
            the molecule's formula.
        molecule (Molecule, optional): The molecule to generate force 
            field parameters for. If not provided, the molecule will be 
            read from the input file.
        mae_cmd (str, optional): The command to run the structconvert 
            utility program. If not provided, the command will be read 
            from the config.ini file.
        ffld_cmd (str, optional): The command to run the ffld_server 
            utility program. If not provided, the command will be read 
            from the config.ini file.
        working_dir (str, optional): The directory to run Maestro in. 
            Defaults to the current working directory.
        maestro_cleanup (bool, optional): Whether to clean up the 
            Maestro files after running. Defaults to False.
        db (str, optional): The connection information of the database 
            to save the force field parameters to. If not provided, the 
            connection information will be read from the db.json file.
        save_ff_to_db (bool, optional): Whether to save the force field 
            parameters to a database. Defaults to False.
        save_ff_to_file (bool, optional): Whether to save the force 
            field parameters to a file. Defaults to True.
        ff_filename (str, optional): The name of the file to save the 
            force field parameters to. Defaults to "ff.json".
        system_force_field_dict (dict, optional): A dictionary to store 
            the force field parameters for multiple molecules. This is 
            only needed if this dictionary is not already in the 
            ``fw_spec``. If this dict is not in the ``fw_spec`` or 
            provided by the user, an empty dict will be created.
    """
    required_params = ["input_file"]
    optional_params = [
        "label",
        "molecule",
        "mae_cmd",
        "ffld_cmd",
        "working_dir",
        "maestro_cleanup",
        "db",
        "save_ff_to_db",
        "save_ff_to_file",
        "ff_filename",
        "system_force_field_dict"
    ]

    def run_task(self, fw_spec):
        input_file = self["input_file"]
        working_dir = self.get("working_dir", os.getcwd())
        db = self.get("db", None)
        ff_filename = self.get("ff_filename", "ff.json")
        save_to_db = self.get("save_ff_to_db", False)
        save_to_file = self.get("save_ff_to_file", True)

        molecule = self.get("molecule", fw_spec["prev_calc_molecule"])

        if not molecule:
            try:
                molecule = process_mol("get_from_file", self["input_file"])
            except Exception as e:
                logger.error("Could not read molecule from file: {}".format(e))
                raise e

        label = self.get("label", get_mol_formula(molecule))

        maestro = MaestroRunner(
            name=label,
            input_file=input_file,
            mae_cmd=self.get("mae_cmd"),
            ffld_cmd=self.get("ffld_cmd"),
            working_dir=working_dir,
        )
        ff_params = maestro.get_opls_params(self.get("maestro_cleanup", False))
        ff_params["Molecule"] = molecule

        if save_to_db:
            ff_doc = ff_params.copy()
            ff_db = get_db(input_db=db)
            ff_db.insert_force_field(ff_doc, "opls", doi=OPLS_DOI)

        if save_to_file:
            ff_doc = process_ff_doc(ff_params.copy(), "opls", doi=OPLS_DOI)
            with open(os.path.join(working_dir, ff_filename), "w") as file:
                json.dump(ff_doc, file)

        labeled_ff_params = add_ff_labels_to_dict(ff_params, label)

        sys_ff_dict = fw_spec.get(
            "system_force_field_dict", self.get("system_force_field_dict", {})
        )
        sys_ff_dict[label] = labeled_ff_params
        return FWAction(update_spec={"system_force_field_dict": sys_ff_dict})


# TODO: tleap Custodian Firetask
