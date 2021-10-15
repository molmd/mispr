# coding: utf-8


# Defines the nuclear magnetic resonance workflow.

import os
import logging

from fireworks import Firework, Workflow

from mispr.gaussian.fireworks.core import CalcFromRunsDBFW
from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.utilities.inputs import handle_gaussian_inputs
from mispr.gaussian.utilities.metadata import get_job_name
from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.gaussian.firetasks.parse_outputs import NMRtoDB

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"


logger = logging.getLogger(__name__)


def get_nmr_tensors(
    mol_operation_type,
    mol,
    db=None,
    name="nmr_tensor_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    nmr_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    **kwargs
):
    fws = []
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)
    gout_keys = ["mol", "mol_nmr"]

    gaussian_inputs = handle_gaussian_inputs(
        {
            "opt": opt_gaussian_inputs,
            "freq": freq_gaussian_inputs,
            "nmr": nmr_gaussian_inputs,
        },
        solvent_gaussian_inputs,
        solvent_properties,
    )
    opt_gaussian_inputs = gaussian_inputs["opt"]
    freq_gaussian_inputs = gaussian_inputs["freq"]
    nmr_gaussian_inputs = gaussian_inputs["nmr"]

    _, label, opt_freq_fws = common_fw(
        mol_operation_type=mol_operation_type,
        mol=mol,
        working_dir=working_dir,
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        gout_key=gout_keys[0],
        skips=skips,
        **kwargs
    )
    fws += opt_freq_fws

    spec = kwargs.pop("spec", {})
    if not skips or len(skips) == 1 and skips[0].lower() == "opt":
        spec.update(
            {"proceed": {"has_gaussian_completed": True, "stationary_type": "Minimum"}}
        )
    else:
        spec.update({"proceed": {"has_gaussian_completed": True}})

    nmr_fw = CalcFromRunsDBFW(
        db,
        input_file="{}_nmr.com".format(label),
        output_file="{}_nmr.out".format(label),
        name=get_job_name(label, "nmr"),
        parents=fws[:],
        gaussian_input_params=nmr_gaussian_inputs,
        working_dir=os.path.join(working_dir, label, "NMR"),
        cart_coords=cart_coords,
        gout_key=gout_keys[1],
        spec=spec,
        **kwargs
    )
    fws.append(nmr_fw)

    fw_analysis = Firework(
        NMRtoDB(
            db=db,
            keys=gout_keys,
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            **{
                i: j
                for i, j in kwargs.items()
                if i in NMRtoDB.required_params + NMRtoDB.optional_params
            }
        ),
        parents=fws[:],
        name="{}-{}".format(label, "nmr_analysis"),
        spec={"_launch_dir": os.path.join(working_dir, label, "analysis")},
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name=get_job_name(label, name),
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS}
    )
