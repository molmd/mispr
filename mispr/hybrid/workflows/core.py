# coding: utf-8


# Defines the DFT-MD hybrid workflow.

import os

from copy import deepcopy

from fireworks import Workflow

from mispr.lammps.workflows.base import lammps_workflow
from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.utilities.inputs import handle_gaussian_inputs
from mispr.gaussian.workflows.base.esp import get_esp_charges
from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Oct 2021"
__version__ = "0.0.1"


def run_hybrid_calcs(
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
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    esp_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    box_data_type="cubic",
    data_file_name="data.mixture",
    analysis_list=None,
    analysis_settings=None,
    **kwargs,
):
    fws = []
    labels = []
    working_dir = working_dir or os.getcwd()
    num_mols = len(mol)
    gout_keys = [f"mol_{i}" for i in range(0, num_mols)]
    if skips is None:
        skips = [None] * num_mols
    if not ff_method:
        ff_method = ["get_from_esp"] * num_mols
        ff_params = [] * num_mols

    mol = recursive_relative_to_absolute_path(mol, working_dir)
    gaussian_inputs = handle_gaussian_inputs(
        {"opt": deepcopy(opt_gaussian_inputs), "freq": deepcopy(freq_gaussian_inputs)},
        solvent_gaussian_inputs,
        solvent_properties,
    )
    opt_gins = gaussian_inputs["opt"]
    freq_gins = gaussian_inputs["freq"]

    process_mol_func = kwargs.pop("process_mol_func", True)
    kwargs.pop("dir_head", None)
    kwargs.pop("dir_structure", None)
    charges = kwargs.pop("charge", [0] * num_mols)

    # either perform geometry optimization and frequency calc for each molecule or
    # the full esp workflow depending on the ff method for each method
    for ind, [operation, molecule, key, skip, molecule_name, charge, ff] in enumerate(
        zip(
            mol_operation_type,
            mol,
            gout_keys,
            skips,
            kwargs.pop("mol_name", [None] * num_mols),
            charges,
            ff_method,
        )
    ):
        if ff == "get_from_esp":
            esp_wf, label = get_esp_charges(
                mol_operation_type=operation,
                mol=molecule,
                db=db,
                working_dir=working_dir,
                opt_gaussian_inputs=deepcopy(opt_gaussian_inputs),
                freq_gaussian_inputs=deepcopy(freq_gaussian_inputs),
                esp_gaussian_inputs=deepcopy(esp_gaussian_inputs),
                solvent_gaussian_inputs=solvent_gaussian_inputs,
                solvent_properties=solvent_properties,
                cart_coords=cart_coords,
                oxidation_states=oxidation_states,
                skips=skip,
                mol_name=molecule_name,
                charge=charge,
                process_mol_func=process_mol_func,
                gout_keys=[key, f"{key}_esp"],
                **kwargs,
            )
            ff_params[ind] = f"{working_dir}/{label}/ESP/{label}_esp"
            fws += esp_wf.fws
        else:
            _, label, opt_freq_init_fws = common_fw(
                mol_operation_type=operation,
                mol=molecule,
                db=db,
                working_dir=working_dir,
                opt_gaussian_inputs=opt_gins,
                freq_gaussian_inputs=freq_gins,
                cart_coords=cart_coords,
                oxidation_states=oxidation_states,
                gout_key=key,
                skips=skip,
                mol_name=molecule_name,
                charge=charge,
                process_mol_func=process_mol_func,
                **kwargs,
            )
            fws += opt_freq_init_fws
        labels.append(label)

    wf = Workflow(
        fws, name=name, **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS}
    )

    # prepare the system_species_data dict in the format required by lammps_workflow
    system_species_data = {}
    for ind, label in enumerate(labels):
        system_species_data[label] = {
            "molecule": gout_keys[ind],
            "molecule_operation_type": "get_from_run_dict",
            "ff_param_method": ff_method[ind],
            "ff_param_data": ff_params[ind],
            "mol_mixture_type": mol_type[ind],
            "mixture_data": mol_data[ind],
            "save_ff_to_db": kwargs.get("save_ff_to_db", False),
            "save_ff_to_file": kwargs.get("save_ff_to_file", True),
            "charge": charges[ind],
        }

    # add the lammps workflow
    md_wf = lammps_workflow(
        system_species_data=system_species_data,
        system_mixture_type=mixture_type,
        box_data=box_data,
        box_data_type=box_data_type,
        data_file_name=data_file_name,
        db=db,
        working_dir=working_dir,
        analysis_list=analysis_list,
        analysis_settings=analysis_settings,
        name=name,
        process_mol_func=False,
        from_fw_spec=True,
        **kwargs,
    )
    wf.append_wf(md_wf, list(wf.id_fw.keys()))
    return wf
