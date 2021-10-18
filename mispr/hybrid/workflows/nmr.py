# coding: utf-8


# Defines the DFT-MD workflow that extracts solvation structures and computes their
# nuclear magnetic resonances.
import os

from copy import deepcopy

from fireworks import Firework, Workflow

from mispr.hybrid.defaults import (
    NMR_GAUSSIAN_INPUTS,
    OPT_GAUSSIAN_INPUTS,
    FREQ_GAUSSIAN_INPUTS,
)
from mispr.hybrid.workflows.core import run_hybrid_calcs
from mispr.hybrid.firetasks.nmr_from_md import NMRFromMD

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Oct 2021"
__version__ = "0.0.1"


def get_solvation_structures_nmr(
    mol_operation_type,
    mol,
    mol_type,
    mol_data,
    box_data,
    ff_method=None,
    ff_params=None,
    mixture_type="number of molecules",
    db=None,
    name="hybrid_calculation",
    working_dir=None,
    opt_esp_gaussian_inputs=None,
    freq_esp_gaussian_inputs=None,
    esp_gaussian_inputs=None,
    opt_nmr_gaussian_inputs=None,
    freq_nmr_gaussian_inputs=None,
    nmr_gaussian_inputs=None,
    esp_solvent_gaussian_inputs=None,
    esp_solvent_properties=None,
    nmr_solvent_gaussian_inputs=None,
    nmr_solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    box_data_type="cubic",
    data_file_name="data.mixture",
    analysis_list=None,
    analysis_settings=None,
    **kwargs,
):
    if not working_dir:
        working_dir = os.getcwd()
    if not analysis_list:
        analysis_list = []
    analysis_list += ["diffusion", "rdf", "cn", "clusters"]
    analysis_list = list(set(analysis_list))

    wf = run_hybrid_calcs(
        mol_operation_type,
        mol,
        mol_type,
        mol_data,
        box_data,
        ff_method,
        ff_params,
        mixture_type,
        db,
        name,
        working_dir,
        opt_esp_gaussian_inputs,
        freq_esp_gaussian_inputs,
        esp_gaussian_inputs,
        esp_solvent_gaussian_inputs,
        esp_solvent_properties,
        cart_coords,
        oxidation_states,
        skips,
        box_data_type,
        data_file_name,
        analysis_list,
        analysis_settings,
        **kwargs,
    )

    dft_md_fw_ids = list(wf.id_fw.keys())

    nmr_dir = f"{working_dir}/nmr"
    nmr_fw = Firework(
        NMRFromMD(
            db=db,
            working_dir=nmr_dir,
            opt_gaussian_inputs=deepcopy(opt_nmr_gaussian_inputs)
            or deepcopy(OPT_GAUSSIAN_INPUTS),
            freq_gaussian_inputs=deepcopy(freq_nmr_gaussian_inputs)
            or deepcopy(FREQ_GAUSSIAN_INPUTS),
            nmr_gaussian_inputs=deepcopy(nmr_gaussian_inputs)
            or deepcopy(NMR_GAUSSIAN_INPUTS),
            solvent_gaussian_inputs=nmr_solvent_gaussian_inputs,
            solvent_properties=nmr_solvent_properties,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            additional_kwargs=kwargs,
        ),
        name="nmr_calculation",
        spec={"_launch_dir": nmr_dir,},
    )
    wf.append_wf(Workflow.from_Firework(nmr_fw), dft_md_fw_ids)
    return wf
