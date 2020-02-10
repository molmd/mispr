import os
import subprocess

from fireworks.core.firework import FiretaskBase, FWAction, Firework
from fireworks.utilities.fw_utilities import explicit_serialize


@explicit_serialize
class RunGaussianDirect(FiretaskBase):
    optional_params = ["working_dir", "input_file", "output_file"]

    def run_task(self, fw_spec):
        working_dir = self.get("working_dir", os.getcwd())

        input_file = self.get('input_file', 'mol.com')
        input_path = os.path.join(working_dir, input_file)

        output_file = self.get('output_file', 'mol.out')
        output_path = os.path.join(working_dir, output_file)

        cmd = f"g09 < {input_path} > {output_path}"
        print("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        print("Finished running with return code: {}".format(return_code))


@explicit_serialize
class RunGaussianCustodian(FiretaskBase):
    def run_task(self, fw_spec):
        pass


@explicit_serialize
class RunGaussianFake(FiretaskBase):
    def run_task(self, fw_spec):
        pass
