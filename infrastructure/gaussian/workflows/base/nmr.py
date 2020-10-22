import os
import logging

from fireworks import Firework, Workflow

from infrastructure.gaussian.utils.utils import get_job_name, \
    get_mol_formula, recursive_relative_to_absolute_path
from infrastructure.gaussian.fireworks.core_standard import CalcFromRunsDBFW
from infrastructure.gaussian.workflows.base.core import common_fw, \
    WORKFLOW_KWARGS

logger = logging.getLogger(__name__)


def get_nmr_tensors(mol_operation_type,
                    mol,
                    db=None,
                    name="nmr_tensor_calculation",
                    working_dir=None,
                    opt_gaussian_inputs=None,
                    freq_gaussian_inputs=None,
                    nmr_gaussian_inputs=None,
                    cart_coords=True,
                    oxidation_states=None,
                    skip_opt_freq=False,
                    **kwargs):
    fws = []
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    _, label, opt_freq_fws = common_fw(mol_operation_type=mol_operation_type,
                                       mol=mol,
                                       working_dir=working_dir,
                                       db=db,
                                       opt_gaussian_inputs=opt_gaussian_inputs,
                                       freq_gaussian_inputs=freq_gaussian_inputs,
                                       cart_coords=cart_coords,
                                       oxidation_states=oxidation_states,
                                       skip_opt_freq=skip_opt_freq,
                                       **kwargs)
    fws += opt_freq_fws
    # mol_formula = get_mol_formula(mol)
    nmr_gaussian_inputs = nmr_gaussian_inputs or {}
    if "route_parameters" not in nmr_gaussian_inputs:
        nmr_gaussian_inputs.update({"route_parameters": {"NMR": "GIAO"}})

    if not skip_opt_freq:
        spec = {"proceed": {"has_gaussian_completed": True,
                            "stationary_type": "Minimum"}}
    else:
        spec = {"proceed": {"has_gaussian_completed": True}}

    nmr_fw = CalcFromRunsDBFW(db,
                              input_file="{}_nmr.com".format(label),
                              output_file="{}_nmr.out".format(label),
                              name=get_job_name(label, "nmr"),
                              parents=fws[:],
                              gaussian_input_params=nmr_gaussian_inputs,
                              working_dir=os.path.join(working_dir, label,
                                                       "NMR"),
                              cart_coords=cart_coords,
                              spec=spec,
                              **kwargs
                              )
    fws.append(nmr_fw)
    return Workflow(fws,
                    name=get_job_name(label, name),
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})
