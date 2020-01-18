import os
import subprocess

from fireworks.core.firework import FiretaskBase, FWAction, Firework
from fireworks.utilities.fw_utilities import explicit_serialize


@explicit_serialize
class RunGaussianDirectTest(FiretaskBase):
    def run_task(self, fw_spec):
        cmd = f"g09 < {fw_spec['input_file_path']} > " \
              f"{os.path.splitext(fw_spec['input_file_path'])[0]}.out"
        print("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        print("Command {} finished running with return code: {}".format(
            cmd, return_code))
        fw_spec['output_path'] = f"{os.path.splitext(fw_spec['input_file_path'])[0]}.out"
        return FWAction(update_spec={'output_path': fw_spec['output_path']})


@explicit_serialize
class RunGaussianCustodian(FiretaskBase):
    def run_task(self, fw_spec):
        pass


@explicit_serialize
class RunGaussianFake(FiretaskBase):
    def run_task(self, fw_spec):
        pass

