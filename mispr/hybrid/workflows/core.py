import os

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.utilities.inputs import handle_gaussian_inputs
from mispr.gaussian.workflows.base.esp import get_esp_charges
from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.lammps.workflows.base import lammps_workflow


def run_hybrid_calcs(
    mol_operation_type,
    mol,
    ff_method,
    ff_params,
    mol_type,
    mol_data,
    box_data,
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
    data_file_name="complex.data",
    analysis_list=None,
    analysis_settings=None,
    **kwargs,
):
    # mol_operation_type = []
    # mol = []
    # system_species_data: [dict]
    #               {
    #               species1_label: {"molecule": pmg.Molecule,
    #                                "ff_param_method": str can be one of the
    #                                     following ("get_from_esp",
    #                                                "get_from_prmtop",
    #                                                "get_from_dict")
    #                                     defaults to "get_from_esp",
    #                                "ff_param_data": str (path_to_file)
    #                                                 or dict (unlabeled_ff_dict),
    #                                "mol_mixture_type": "Solutes" or "Solvents",
    #                                "mixture_data": int or dict depends on
    #                                     system_mixture_type},
    #               ...,
    #               speciesN_label: {...}
    #               }

    fws = []
    labels = []
    working_dir = working_dir or os.getcwd()
    num_mols = len(mol)
    gout_keys = [f"mol_{i}" for i in range(1, num_mols + 1)]
    if skips is None:
        skips = [None] * num_mols

    mol = recursive_relative_to_absolute_path(mol, working_dir)
    gaussian_inputs = handle_gaussian_inputs(
        {"opt": opt_gaussian_inputs, "freq": freq_gaussian_inputs},
        solvent_gaussian_inputs,
        solvent_properties,
    )
    opt_gaussian_inputs = gaussian_inputs["opt"]
    freq_gaussian_inputs = gaussian_inputs["freq"]

    parents = []

    # perform geometry optimization and frequency calc for each molecule
    for [operation, molecule, key, skip, molecule_name, charge] in zip(
        mol_operation_type,
        mol,
        gout_keys,
        skips,
        kwargs.pop("mol_name", [None] * num_mols),
        kwargs.pop("charge", [0] * num_mols),
    ):
        _, label, opt_freq_init_fws = common_fw(
            mol_operation_type=operation,
            mol=molecule,
            working_dir=working_dir,
            db=db,
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            gout_key=key,
            skips=skip,
            mol_name=molecule_name,
            charge=charge,
            **kwargs,
        )
        fws += opt_freq_init_fws
        parents.append(len(fws))
        labels.append(label)

    wf = Workflow(
        fws,
        name=name,
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS}
    )

    # perform esp calc for molecules that user requests amber ff for
    # TODO: check if all kwargs can be given to this workflow
    # TODO: correctly add the parent-child relationship; currently doing it after all opt_freq calcs
    for ind, ff in enumerate(ff_method):
        if ff == "get_from_esp":
            esp_wf = get_esp_charges(
                mol_operation_type="get_from_run_dict",
                mol=gout_keys[ind + 1],
                db=db,
                working_dir=working_dir,
                opt_gaussian_inputs=opt_gaussian_inputs,
                freq_gaussian_inputs=freq_gaussian_inputs,
                esp_gaussian_inputs=esp_gaussian_inputs,
                solvent_gaussian_inputs=solvent_gaussian_inputs,
                solvent_properties=solvent_properties,
                cart_coords=cart_coords,
                oxidation_states=oxidation_states,
                skips=["opt", "freq"],
                mol_name=labels[ind + 1],
                **kwargs,
            )
            wf.append_wf(esp_wf, fws[:])

    # firetask/function that prepares inputs to the lammps workflow
    # system_species_data = {}
    # md_wf = lammps_workflow(system_species_data=system_species_data,
    #                         system_mixture_type=mol_type,
    #                         box_data=box_data, box_data_type=box_data_type,
    #                         db=db,
    #                         working_dir=working_dir,
    #                         analysis_list=analysis_list,
    #                         analysis_settings=analysis_settings,
    #                         name=name,
    #                         **kwargs)
    # wf.append(md_wf, wf.fws[-1])
    return wf
