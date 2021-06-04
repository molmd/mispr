# coding: utf-8


# Defines firetasks for running Gaussian calculations.

import os
import logging
import subprocess

from configparser import ConfigParser

from timeit import default_timer as timer

from fireworks.core.firework import FiretaskBase
from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks.utilities.fw_utilities import explicit_serialize

from custodian import Custodian
from custodian.gaussian.jobs import GaussianJob
from custodian.gaussian.handlers import GaussianErrorHandler

from infrastructure.gaussian.defaults import CUSTODIAN_MAX_ERRORS

__author__ = 'Rasha Atwi'
__maintainer__ = 'Rasha Atwi'
__email__ = 'rasha.atwi@stonybrook.edu'
__status__ = 'Development'
__date__ = 'Jan 2021'
__version__ = 0.2

logger = logging.getLogger(__name__)


@explicit_serialize
class RunGaussianDirect(FiretaskBase):
    required_params = []
    optional_params = ['input_file', 'output_file', 'gaussian_cmd']

    def run_task(self, fw_spec):
        working_dir = os.getcwd()

        input_file = self.get('input_file', 'mol.com')
        input_path = os.path.join(working_dir, input_file)

        output_file = self.get('output_file', 'mol.out')
        output_path = os.path.join(working_dir, output_file)

        cmd = self.get('gaussian_cmd')
        if not cmd:
            cfg = ConfigParser()
            cfg.read(CONFIG_FILE_DIR + '/config.ini')
            cmd = cfg['RunCalc']['gcmd']
        cmd = cmd.replace('$input_path$', input_path). \
            replace('$output_path$', output_path)

        logger.info('Running command: {}'.format(cmd))
        st = timer()
        return_code = subprocess.call(cmd, shell=True)
        run_time = timer() - st
        logger.info('Finished running with return code: {}'.format(return_code))
        fw_spec['run_time'] = run_time


@explicit_serialize
class RunGaussianCustodian(FiretaskBase):
    required_params = []
    optional_params = ['input_file', 'output_file', 'gaussian_cmd',
                       'stderr_file', 'job_type', 'backup', 'scf_max_cycles',
                       'opt_max_cycles', 'cart_coords', 'max_errors',
                       'scf_functional', 'scf_basis_set', 'prefix',
                       'check_convergence']

    def run_task(self, fw_spec):
        # working_dir = os.getcwd()

        input_file = self.get('input_file', 'mol.com')
        # input_path = os.path.join(working_dir, input_file)

        output_file = self.get('output_file', 'mol.out')
        # output_path = os.path.join(working_dir, output_file)

        backup = self.get('backup', True)
        stderr_file = self.get('stderr_file', 'stderr.txt')
        scf_max_cycles = self.get('scf_max_cycles', 100)
        opt_max_cycles = self.get('opt_max_cycles', 100)
        max_errors = self.get('max_errors', CUSTODIAN_MAX_ERRORS)

        job_type = self.get('job_type', 'normal')
        scf_functional = self.get('scf_functional', None)
        scf_basis_set = self.get('scf_basis_set', None)
        cart_coords = self.get('cart_coords', True)
        check_convergence = self.get('check_convergence', True)

        cmd = self.get('gaussian_cmd')
        if not cmd:
            cfg = ConfigParser()
            cfg.read(CONFIG_FILE_DIR + '/config.ini')
            cmd = cfg['RunCalc']['gcmd']
        cmd = cmd.replace('$input_path$', input_file). \
            replace('$output_path$', output_file)

        if job_type == 'normal':
            jobs = [GaussianJob(gaussian_cmd=cmd,
                                input_file=input_file,
                                output_file=output_file,
                                stderr_file=stderr_file,
                                backup=backup)]

        elif job_type == 'better_scf_guess':
            if not scf_functional or not scf_basis_set:
                raise Exception(f'{job_type} is requested but the functional '
                                f'and/or basis set to use for the SCF '
                                f'calculation are not provided! Exiting...')
            jobs = GaussianJob.better_scf_guess(gaussian_cmd=cmd,
                                                input_file=input_file,
                                                output_file=output_file,
                                                stderr_file=stderr_file,
                                                backup=backup,
                                                cart_coords=cart_coords)
        else:
            raise ValueError(f'Unsupported job type: {job_type}')

        handlers = [GaussianErrorHandler(input_file=input_file,
                                         output_file=output_file,
                                         stderr_file=stderr_file,
                                         cart_coords=cart_coords,
                                         scf_max_cycles=scf_max_cycles,
                                         opt_max_cycles=opt_max_cycles,
                                         job_type=job_type,
                                         scf_functional=scf_functional,
                                         scf_basis_set=scf_basis_set,
                                         prefix=self.get('prefix', 'error'),
                                         check_convergence=check_convergence)]

        c = Custodian(handlers,
                      jobs,
                      max_errors=max_errors)
        st = timer()
        return_code = c.run()
        run_time = timer() - st
        logger.info('Finished running with return code: {}'.format(return_code))
        fw_spec['run_time'] = run_time


@explicit_serialize
class RunGaussianFake(FiretaskBase):
    def run_task(self, fw_spec):
        pass
