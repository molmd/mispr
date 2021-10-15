# coding: utf-8


# Defines a list of common fireworks used in Gaussian workflows.

import os
import logging

from copy import deepcopy
from fireworks import Firework, Workflow

from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.fireworks.core import CalcFromMolFW, CalcFromRunsDBFW
from mispr.gaussian.utilities.gout import process_run
from mispr.gaussian.utilities.metadata import get_job_name, get_mol_formula
from mispr.gaussian.firetasks.parse_outputs import ProcessRun

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

WORKFLOW_KWARGS = Workflow.__init__.__code__.co_varnames


class GoutTypeError(Exception):
    def __init__(self, msg=None):
        if msg is None:
            msg = (
                "no Gaussian output provided; to skip optimization, you "
                "need to input the molecule in any of the "
                "supported Gaussian output formats"
            )
        self.msg = msg
        super(GoutTypeError, self).__init__(msg)

    def __repr__(self):
        return self.msg


def _recursive_key_check(grun, grun_key):
    key_exists = grun_key in grun
    if not key_exists:
        for key, value in grun.items():
            if isinstance(value, dict):
                key_exists = key_exists or _recursive_key_check(value, grun_key)
    return key_exists


def _process_mol_check(
    working_dir,
    process_mol_func=True,
    mol_operation_type=None,
    mol=None,
    dir_head=None,
    mol_name=None,
    db=None,
    dir_structure=None,
    charge=None,
):
    if process_mol_func:
        mol = process_mol(mol_operation_type, mol, db=db, charge=charge,
                          working_dir=working_dir)
        mol_operation_type = "get_from_mol"
        mol_formula = get_mol_formula(mol)
        opt_job_name = get_job_name(mol, "optimization")
        freq_job_name = get_job_name(mol, "frequency")
        label = mol_formula
    elif mol_name:
        opt_job_name = "{}_optimization".format(mol_name)
        freq_job_name = "{}_frequency".format(mol_name)
        label = mol_name
    else:
        opt_job_name = "optimization"
        freq_job_name = "frequency"
        label = "mol"
    if not dir_head:
        dir_head = label
    # dir_struct = [dir_head] + kwargs.get('dir_structure', [])
    if dir_structure is None:
        dir_structure = []
    dir_struct = [dir_head] + dir_structure
    working_dir = os.path.join(working_dir, *dir_struct)
    return mol_operation_type, mol, label, opt_job_name, freq_job_name, working_dir


# TODO: avoid overwriting a directory if it exists (happens when molecules have
#  same mol formula)
def common_fw(
    mol_operation_type,
    mol,
    working_dir,
    opt_gaussian_inputs,
    freq_gaussian_inputs,
    cart_coords,
    oxidation_states,
    gout_key=None,
    db=None,
    process_mol_func=True,
    mol_name=None,
    dir_head=None,
    skips=None,
    check_result=None,
    **kwargs,
):
    fws = []
    if not gout_key:
        gout_key = "mol"

    original_mol_operation_type = deepcopy(mol_operation_type)
    original_mol = deepcopy(mol)

    (
        mol_operation_type,
        mol,
        label,
        opt_job_name,
        freq_job_name,
        working_dir,
    ) = _process_mol_check(
        working_dir,
        process_mol_func,
        mol_operation_type,
        mol,
        dir_head,
        mol_name,
        db,
        kwargs.get("dir_structure", []),
        kwargs.get("charge", None)
    )

    if skips is None:
        skips = []
    if not skips:
        # if user chooses to perform both opt and freq calc
        opt_fw = CalcFromMolFW(
            mol=mol,
            mol_operation_type=mol_operation_type,
            db=db,
            name=opt_job_name,
            working_dir=os.path.join(working_dir, "Optimization"),
            input_file=f"{label}_opt.com",
            output_file=f"{label}_opt.out",
            gaussian_input_params=opt_gaussian_inputs,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            gout_key=gout_key + "_opt",
            **kwargs,
        )
        fws.append(opt_fw)

        spec = kwargs.pop("spec", {})
        spec.update({"proceed": {"has_gaussian_completed": True}})
        freq_fw = CalcFromRunsDBFW(
            db=db,
            name=freq_job_name,
            parents=opt_fw,
            gaussian_input_params=freq_gaussian_inputs,
            working_dir=os.path.join(working_dir, "Frequency"),
            input_file=f"{label}_freq.com",
            output_file=f"{label}_freq.out",
            cart_coords=cart_coords,
            gout_key=gout_key,
            spec=spec,
            **kwargs,
        )
        fws.append(freq_fw)

    elif len(skips) == 1:
        if skips[0].lower() == "opt":
            # if user chooses to skip opt, only perform freq calc and restrict
            # the mol input to any gaussian output format
            if mol_operation_type == "get_from_mol" and len(mol) == 1:
                freq_fw = CalcFromMolFW(
                    mol=mol,
                    mol_operation_type=mol_operation_type,
                    db=db,
                    name=freq_job_name,
                    working_dir=os.path.join(working_dir, "Frequency"),
                    input_file=f"{label}_freq.com",
                    output_file=f"{label}_freq.out",
                    gaussian_input_params=freq_gaussian_inputs,
                    cart_coords=cart_coords,
                    oxidation_states=oxidation_states,
                    gout_key=gout_key,
                    **kwargs,
                )
                fws.append(freq_fw)
            elif original_mol_operation_type not in [
                "get_from_gout",
                "get_from_gout_file",
                "get_from_run_dict",
                "get_from_run_id",
                "get_from_run_query",
            ]:
                raise GoutTypeError()
            else:
                freq_fw = CalcFromMolFW(
                    mol=mol,
                    mol_operation_type=mol_operation_type,
                    db=db,
                    name=freq_job_name,
                    working_dir=os.path.join(working_dir, "Frequency"),
                    input_file=f"{label}_freq.com",
                    output_file=f"{label}_freq.out",
                    gaussian_input_params=freq_gaussian_inputs,
                    cart_coords=cart_coords,
                    oxidation_states=oxidation_states,
                    gout_key=gout_key,
                    **kwargs,
                )
                fws.append(freq_fw)

        elif skips[0].lower() == "freq":
            # if user chooses to skip freq, only perform opt
            opt_fw = CalcFromMolFW(
                mol=mol,
                mol_operation_type=mol_operation_type,
                db=db,
                name=opt_job_name,
                working_dir=os.path.join(working_dir, "Optimization"),
                input_file=f"{label}_opt.com",
                output_file=f"{label}_opt.out",
                gaussian_input_params=opt_gaussian_inputs,
                cart_coords=cart_coords,
                oxidation_states=oxidation_states,
                gout_key=gout_key + "_opt",
                **kwargs,
            )
            fws.append(opt_fw)

    else:
        # if user chooses to skip both opt and freq, only process the mol and
        # restrict the mol input to any gaussian output format
        if original_mol_operation_type not in [
            "get_from_gout",
            "get_from_gout_file",
            "get_from_run_dict",
            "get_from_run_id",
            "get_from_run_query",
        ]:
            raise GoutTypeError()
        else:
            run = process_run(original_mol_operation_type, original_mol, db=db,
                              working_dir=working_dir)
            if check_result:
                keys_exist = []
                for grun_key in check_result:
                    keys_exist.append(_recursive_key_check(run, grun_key))
                not_found_ind = [
                    ind for ind, boolean in enumerate(keys_exist) if not boolean
                ]
                if not_found_ind:
                    not_found_keys = [check_result[ind] for ind in not_found_ind]
                    raise ValueError(
                        "Gaussian output does not include {}. "
                        "Stopping.".format(not_found_keys)
                    )

            spec = kwargs.pop("spec", {})
            spec.update({"tag": kwargs.get("tag", "unknown")})
            spec.update({"_launch_dir": working_dir})
            fws.append(
                Firework(
                    ProcessRun(
                        run=run,
                        operation_type="get_from_run_dict",
                        db=db,
                        gout_key=gout_key,
                        **{
                            i: j
                            for i, j in kwargs.items()
                            if i
                            in ProcessRun.required_params + ProcessRun.optional_params
                        },
                    ),
                    name=get_job_name(label, "process_run"),
                    spec=spec,
                )
            )
    return mol, label, fws
