# coding: utf-8


# Defines firetasks for running Gaussian calculations.

import os
import shutil
import logging
import subprocess

from timeit import default_timer as timer
from configparser import ConfigParser

import numpy as np

from monty.os.path import zpath
from monty.serialization import loadfn

from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks.core.firework import Firework, FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from pymatgen.io.gaussian import GaussianInput

from mispr.gaussian.utilities.misc import recursive_compare_dicts

from custodian import Custodian
from custodian.gaussian.jobs import GaussianJob
from custodian.gaussian.handlers import WalTimeErrorHandler, GaussianErrorHandler

from mispr.gaussian.defaults import CUSTODIAN_MAX_ERRORS

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)


@explicit_serialize
class RunGaussianDirect(FiretaskBase):
    required_params = []
    optional_params = ["input_file", "output_file", "gaussian_cmd"]

    def run_task(self, fw_spec):
        working_dir = os.getcwd()

        input_file = self.get("input_file", "mol.com")
        input_path = os.path.join(working_dir, input_file)

        output_file = self.get("output_file", "mol.out")
        output_path = os.path.join(working_dir, output_file)

        cmd = self.get("gaussian_cmd")
        if not cmd:
            cfg = ConfigParser()
            cfg.read(CONFIG_FILE_DIR + "/config.ini")
            cmd = cfg["RunCalc"]["gcmd"]
        cmd = cmd.replace("$input_path$", input_path).replace(
            "$output_path$", output_path
        )

        logger.info("Running command: {}".format(cmd))
        st = timer()
        return_code = subprocess.call(cmd, shell=True)
        run_time = timer() - st
        logger.info("Finished running with return code: {}".format(return_code))
        fw_spec["run_time"] = run_time


@explicit_serialize
class RunGaussianCustodian(FiretaskBase):
    required_params = []
    optional_params = [
        "input_file",
        "output_file",
        "gaussian_cmd",
        "stderr_file",
        "job_type",
        "backup",
        "scf_max_cycles",
        "opt_max_cycles",
        "cart_coords",
        "max_errors",
        "lower_functional",
        "lower_basis_set",
        "prefix",
        "suffix",
        "check_convergence",
        "wall_time",
        "buffer_time",
        "max_wall_time_corrections",
        "additional_fw",
    ]

    def run_task(self, fw_spec):
        # working_dir = os.getcwd()

        input_file = self.get("input_file", "mol.com")
        wt_input_file = fw_spec.get("overwrite_input_file", input_file)

        # input_path = os.path.join(working_dir, input_file)

        output_file = self.get("output_file", "mol.out")
        # output_path = os.path.join(working_dir, output_file)

        backup = self.get("backup", True)
        prefix = self.get("prefix", "error")
        stderr_file = self.get("stderr_file", "stderr.txt")
        scf_max_cycles = self.get("scf_max_cycles", 100)
        opt_max_cycles = self.get("opt_max_cycles", 100)
        max_errors = self.get("max_errors", CUSTODIAN_MAX_ERRORS)

        job_type = self.get("job_type", "normal")
        lower_functional = self.get("lower_functional", None)
        lower_basis_set = self.get("lower_basis_set", None)
        cart_coords = self.get("cart_coords", True)
        check_convergence = self.get("check_convergence", True)
        wall_time = self.get("wall_time", None)
        buffer_time = self.get("buffer_time", 300)

        cmd = self.get("gaussian_cmd")
        if not cmd:
            cfg = ConfigParser()
            cfg.read(CONFIG_FILE_DIR + "/config.ini")
            cmd = cfg["RunCalc"]["gcmd"]
        cmd = cmd.replace("$input_path$", wt_input_file).replace(
            "$output_path$", output_file
        )

        if job_type == "normal":
            jobs = [
                GaussianJob(
                    gaussian_cmd=cmd,
                    input_file=input_file,
                    output_file=output_file,
                    stderr_file=stderr_file,
                    suffix=self.get("suffix", ""),
                    backup=backup,
                )
            ]

        elif job_type == "better_guess":
            if not lower_functional or not lower_basis_set:
                raise Exception(
                    f"{job_type} is requested but the functional "
                    f"and/or basis set to use for the SCF "
                    f"calculation are not provided! Exiting..."
                )
            jobs = GaussianJob.better_guess(
                gaussian_cmd=cmd,
                input_file=input_file,
                output_file=output_file,
                stderr_file=stderr_file,
                backup=backup,
                cart_coords=cart_coords,
            )
        else:
            raise ValueError(f"Unsupported job type: {job_type}")

        handlers = [
            GaussianErrorHandler(
                input_file=input_file,
                output_file=output_file,
                stderr_file=stderr_file,
                cart_coords=cart_coords,
                scf_max_cycles=scf_max_cycles,
                opt_max_cycles=opt_max_cycles,
                job_type=job_type,
                lower_functional=lower_functional,
                lower_basis_set=lower_basis_set,
                prefix=prefix,
                check_convergence=check_convergence,
            )
        ]
        if wall_time:
            handlers.append(
                WalTimeErrorHandler(
                    wall_time=wall_time,
                    buffer_time=buffer_time,
                    input_file=input_file,
                    output_file=output_file,
                    stderr_file=stderr_file,
                    prefix=prefix,
                )
            )

        c = Custodian(handlers, jobs, max_errors=max_errors)

        st = timer()
        try:
            return_code = c.run()
        except Exception as e:
            # TODO: add a checkpoint here that fw_id is accessible
            #  (only if _add_launchpad_and_fw_id is True in the fw_spec)
            if (
                os.path.exists(zpath("custodian.json"))
                and os.path.getsize("custodian.json") > 0
            ):
                custodian_data = loadfn(zpath("custodian.json"))
                for entry in custodian_data:
                    for correction in entry.get("corrections"):
                        if "wall_time_limit" in correction.get("errors"):
                            if fw_spec.get(
                                "number_of_wall_time_corrections", 0
                            ) <= self.get("max_wall_time_corrections", 3):
                                print(
                                    "correction number:",
                                    fw_spec.get("number_of_wall_time_corrections", 0),
                                )
                                fw = self.launchpad.get_fw_by_id(self.fw_id)
                                fw.spec.update(
                                    {
                                        "_recovery": self.launchpad.get_recovery(
                                            self.fw_id
                                        ),
                                        "number_of_wall_time_corrections": fw_spec.get(
                                            "number_of_wall_time_corrections", 0
                                        )
                                        + 1,
                                        "overwrite_input_file": input_file + ".wt",
                                    }
                                )
                                new_fw = Firework(fw.tasks, fw.spec, fw.name)
                                return FWAction(detours=[new_fw])
            raise e
        run_time = timer() - st
        logger.info("Finished running with return code: {}".format(return_code))
        fw_spec["run_time"] = run_time


@explicit_serialize
class RunGaussianFake(FiretaskBase):
    required_params = ["ref_dir"]
    optional_params = ["working_dir", "input_file", "tolerance"]

    def run_task(self, fw_spec):
        self._verify_inputs()
        self._clear_inputs()
        self._generate_outputs()

    @staticmethod
    def _recursive_lowercase(obj):
        if isinstance(obj, dict):
            updated_obj = {}
            for k, v in obj.items():
                updated_obj[k.lower()] = RunGaussianFake._recursive_lowercase(v)
            return updated_obj
        elif isinstance(obj, str):
            return obj.lower()
        elif hasattr(obj, "__iter__"):
            updated_obj = []
            for i in obj:
                updated_obj.append(RunGaussianFake._recursive_lowercase(i))
            return updated_obj
        else:
            return obj

    def _verify_inputs(self):
        ref_dir = self["ref_dir"]
        working_dir = self.get("working_dir", os.getcwd())
        gin_file = self.get("input_file", "mol.com")

        user_gin = GaussianInput.from_file(f"{working_dir}/{gin_file}")
        ref_gin = GaussianInput.from_file(f"{ref_dir}/{gin_file}")
        tol = self.get("tolerance") or 0.0001

        np.testing.assert_equal(ref_gin.molecule.species, user_gin.molecule.species)
        np.testing.assert_allclose(
            ref_gin.molecule.cart_coords, user_gin.molecule.cart_coords, atol=tol
        )

        ref_dict = self._recursive_lowercase(ref_gin.as_dict())
        del ref_dict["molecule"]
        user_dict = self._recursive_lowercase(user_gin.as_dict())
        del user_dict["molecule"]
        diff = recursive_compare_dicts(ref_dict, user_dict, "ref_dict", "user_dict")

        if diff:
            raise ValueError(
                f"Gaussian input is inconsistent with reference input!\n{diff}!"
            )
        logger.info("RunGausianFake: verified input successfully")

    def _clear_inputs(self):
        working_dir = self.get("working_dir", os.getcwd())
        gin_file = self.get("input_file", "mol.com")
        gin_file = f"{working_dir}/{gin_file}"
        if os.path.exists(gin_file):
            os.remove(gin_file)

    def _generate_outputs(self):
        ref_dir = self["ref_dir"]
        working_dir = self.get("working_dir", os.getcwd())
        for file in os.listdir(ref_dir):
            full_path = f"{ref_dir}/{file}"
            if os.path.isfile(full_path):
                shutil.copy(full_path, working_dir)
        logger.info("RunGaussianFake: ran fake Gaussian, generated outputs")
