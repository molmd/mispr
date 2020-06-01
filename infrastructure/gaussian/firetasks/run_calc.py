import os
import logging
import subprocess

from configparser import ConfigParser
from timeit import default_timer as timer

from fireworks.core.firework import FiretaskBase
from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks.utilities.fw_utilities import explicit_serialize

logger = logging.getLogger(__name__)


@explicit_serialize
class RunGaussianDirect(FiretaskBase):
    required_params = []
    optional_params = ["input_file", "output_file", "gaussian_cmd"]

    def run_task(self, fw_spec):
        working_dir = os.getcwd()

        input_file = self.get('input_file', 'mol.com')
        input_path = os.path.join(working_dir, input_file)

        output_file = self.get('output_file', 'mol.out')
        output_path = os.path.join(working_dir, output_file)

        cmd = self.get("gaussian_cmd")
        if not cmd:
            cfg = ConfigParser()
            cfg.read(CONFIG_FILE_DIR+'/config.ini')
            cmd = cfg['RunCalc']['gcmd']
        cmd = cmd.replace('$input_path$', input_path).\
            replace('$output_path$', output_path)

        logger.info("Running command: {}".format(cmd))
        st = timer()
        return_code = subprocess.call(cmd, shell=True)
        run_time = timer() - st
        logger.info("Finished running with return code: {}".format(return_code))
        fw_spec["run_time"] = run_time


@explicit_serialize
class RunGaussianCustodian(FiretaskBase):
    def run_task(self, fw_spec):
        pass


@explicit_serialize
class RunGaussianFake(FiretaskBase):
    def run_task(self, fw_spec):
        pass
