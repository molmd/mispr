# coding: utf-8


# Defines firetasks for parsing Gaussian output files.
import os
import logging
import json

from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize


@explicit_serialize
class GaussianToDB(FiretaskBase):
    """
    Assumes the functional/basis_set to be the first parameters on the
    route_parameters line.
    """
    optional_params = [
        "working_dir", "calc_loc", "input_file", "output_file",
        "additional_fields", "db_file", "fw_spec_field", "multirun"
    ]
    # required_params = []
    # optional_params = ["output_path"]

    def run_task(self, fw_spec):
        from pymatgen.io.gaussian import GaussianOutput, GaussianInput
        from pymatgen.core.structure import Molecule
        from infrastructure.gaussian.database import GaussianCalcDb

        working_dir = \
            fw_spec.get('working_dir', self.get("working_dir", os.getcwd()))

        input_file = fw_spec.get('input_file',
                                 self.get("input_file", "mol.com"))
        output_file = fw_spec.get('output_file',
                                  self.get("output_file", "mol.out"))
        #TODO: Include multirun format
        #multirun = self.get("multirun", False)

        logger.info("PARSING DIRECTORY: {}".format(working_dir))

        output_path = os.path.join(working_dir, output_file)
        #TODO: modify to allow taking an output file only (in cases we do not
        # have the input
        #TODO: modify main pymatgen function to take link0_parameters and dieze_tag
        # from the output file
        gout = GaussianOutput(output_path).as_dict()
        # if 'Mulliken_charges':
        #     gout['Mulliken_charges'] = \
        #         json.loads(json.dumps(gout['Mulliken_charges']))
        # if "IOp(^(\d?[1-9]|[1-9]0)$ / " \
        #    "[99-181])":
        #     gout['input']['route'] =

        # Save input file to database in a standard way
        # Save output file to database in
        # Easy searching for a molecule or a document and passing to a firetask,
        # firework, and workflow


