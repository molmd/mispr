# coding: utf-8


# Defines firetasks for parsing Gaussian output files.
import os
import logging
import json

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize

logger = logging.getLogger(__name__)

JOB_TYPES = {'SP', 'Opt', 'Freq', 'IRC', 'IRCMax', 'Scan', 'Polar', 'ADMP',
             'BOMD', 'EET', 'Force', 'Stable', 'Volume', 'Density', 'Guess',
             'Pop', 'SCRF', 'CPHF', 'Prop', 'NMR', 'CIS', 'ZIndo', 'TD', 'EOM',
             'SAC-CI'}


@explicit_serialize
class GaussianToDB(FiretaskBase):
    """
    Assumes the functional/basis_set to be the first parameters on the
    route_parameters line.
    """
    required_params = ["db"]
    optional_params = [
        "working_dir", "input_file", "output_file", "fw_spec_field"
    ]

    def run_task(self, fw_spec):
        from pymatgen.io.gaussian import GaussianOutput, GaussianInput
        from pymatgen.core.structure import Molecule
        from infrastructure.gaussian.database import GaussianCalcDb
        # TODO: Include multirun format
        # TODO: modify to allow taking an output file only (in cases we do not
        # have the input
        working_dir = self.get("working_dir", os.getcwd())
        input_file = self.get("input_file", "mol.com")
        output_file = self.get("output_file", "mol.out")

        logger.info("Parsing directory: {}".format(working_dir))

        output_path = os.path.join(working_dir, output_file)
        gout = GaussianOutput(output_path).as_dict()

        input_path = os.path.join(working_dir, input_file)
        gin = GaussianInput.from_file(input_path).as_dict()

        gout['input']['input_parameters'] = gin['input_parameters']
        task_doc = {'input': gin, 'output': gout}

        if "tag" in fw_spec:
            task_doc.update({"tag": fw_spec["tag"]})

        # Check for additional keys to set based on the fw_spec
        if self.get("fw_spec_field"):
            task_doc.update({self.get("fw_spec_field"):
                                 fw_spec.get(self.get("fw_spec_field"))})

        task_doc['smiles'] = \
            GaussianCalcDb.get_smiles(Molecule.from_dict(gin['molecule']))
        task_doc['functional'] = gin['functional']
        task_doc['basis'] = gin['basis_set']

        job_types = \
            list(filter(lambda x: x in gin['route_parameters'], JOB_TYPES))
        task_doc['type'] = ';'.join(job_types)
        task_doc = json.loads(json.dumps(task_doc))
        print(task_doc)

        runs_db = GaussianCalcDb(**fw_spec['db'])
        runs_db.insert_run(task_doc)
        logger.info("Finished parsing output and saving to db")
        return FWAction(update_spec={"gaussian_output": task_doc})


