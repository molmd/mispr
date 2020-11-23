import os
import logging

from fireworks import Firework, Workflow

from infrastructure.gaussian.utils.utils import get_job_name, \
    recursive_relative_to_absolute_path
from infrastructure.gaussian.fireworks.core_standard import CalcFromRunsDBFW
from infrastructure.gaussian.workflows.base.core import common_fw, \
    WORKFLOW_KWARGS

logger = logging.getLogger(__name__)


def get_esp_charges(mol_operation_type,
                    mol,
                    db=None,
                    name="esp_charges_calculation",
                    working_dir=None,
                    opt_gaussian_inputs=None,
                    freq_gaussian_inputs=None,
                    esp_gaussian_inputs=None,
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

    esp_gaussian_inputs = esp_gaussian_inputs or {}
    if "route_parameters" not in esp_gaussian_inputs:
        esp_gaussian_inputs.update({"route_parameters": {"pop": "MK",
                                                         "iop(6/50=1)": None}})
    # input_parameters from a previous run are overwritten
    if "input_parameters" not in esp_gaussian_inputs:
        mol_esp = os.path.join(
            working_dir, "{}_esp".format(
                os.path.join(working_dir, label, "ESP", label)))
        esp_gaussian_inputs.update({"input_parameters": {mol_esp: None}})

    if not skip_opt_freq:
        spec = {"proceed": {"has_gaussian_completed": True,
                            "stationary_type": "Minimum"}}
    else:
        spec = {"proceed": {"has_gaussian_completed": True}}

    esp_fw = CalcFromRunsDBFW(db,
                              input_file="{}_esp.com".format(label),
                              output_file="{}_esp.out".format(label),
                              name=get_job_name(label, "esp"),
                              parents=fws[:],
                              gaussian_input_params=esp_gaussian_inputs,
                              working_dir=os.path.join(working_dir, label,
                                                       "ESP"),
                              cart_coords=cart_coords,
                              spec=spec,
                              **kwargs
                              )
    fws.append(esp_fw)
    return Workflow(fws,
                    name=get_job_name(label, name),
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})
