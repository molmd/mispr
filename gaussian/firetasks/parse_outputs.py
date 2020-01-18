# coding: utf-8


# Defines firetasks for parsing Gaussian output files.

from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize


@explicit_serialize
class GaussianToDB(FiretaskBase):
    required_params = []
    optional_params = ["output_path"]

    def run_task(self, fw_spec):
        # TODO:
        from pymatgen.io.gaussian import GaussianOutput
        out = GaussianOutput(fw_spec["output_path"]).as_dict()

        # Save input file to database in a standard way
        # Save output file to database in
        # Easy searching for a molecule or a document and passing to a firetask,
        # firework, and workflow


