# coding: utf-8


# Defines the bond dissociation energy workflow.

import os

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.utilities.inputs import handle_gaussian_inputs
from mispr.gaussian.utilities.metadata import get_job_name
from mispr.gaussian.fireworks.break_mol import BreakMolFW
from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.gaussian.firetasks.parse_outputs import BDEtoDB

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"


def get_bde(
    mol_operation_type,
    mol,
    ref_charge=0,
    fragment_charges=None,
    bonds=None,
    open_rings=False,
    db=None,
    name="bde_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=False,
    visualize=True,
    **kwargs
):
    fws = []
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)
    gout_key = "ref_mol"
    gaussian_inputs = handle_gaussian_inputs(
        {"opt": opt_gaussian_inputs, "freq": freq_gaussian_inputs},
        solvent_gaussian_inputs,
        solvent_properties,
    )
    opt_gaussian_inputs = gaussian_inputs["opt"]
    freq_gaussian_inputs = gaussian_inputs["freq"]

    if skips:
        check_result = ["final_energy", "Enthalpy"]
    else:
        check_result = None

    _, label, opt_freq_fws = common_fw(
        mol_operation_type=mol_operation_type,
        mol=mol,
        charge=ref_charge,
        working_dir=working_dir,
        dir_structure=["principle_mol"],
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        skips=skips,
        check_result=check_result,
        gout_key=gout_key,
        **kwargs
    )
    fws += opt_freq_fws

    break_fw = BreakMolFW(
        mol=gout_key,
        mol_operation_type="get_from_run_dict",
        from_fw_spec=True,
        bonds=bonds,
        open_rings=open_rings,
        ref_charge=ref_charge,
        fragment_charges=fragment_charges,
        db=db,
        calc_frags=True,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        cart_coords=cart_coords,
        name=get_job_name(label, "breaking"),
        parents=fws[:],
        working_dir=os.path.join(working_dir, label, "fragments"),
        **kwargs
    )
    fws.append(break_fw)

    fw_analysis = Firework(
        BDEtoDB(
            principle_mol_key=gout_key,
            db=db,
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            visualize=visualize,
            **{
                i: j
                for i, j in kwargs.items()
                if i in BDEtoDB.required_params + BDEtoDB.optional_params
            }
        ),
        parents=fws[:],
        name="{}-{}".format(label, "bde_analysis"),
        spec={
            "_launch_dir": os.path.join(working_dir, label, "analysis"),
            "_allow_fizzled_parents": True,
        },
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name="{}_{}".format(label, name),
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS}
    )
