# coding: utf-8


# Defines the binding energy workflow.

import os
import logging

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.utilities.inputs import handle_gaussian_inputs
from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.gaussian.firetasks.parse_outputs import BindingEnergytoDB

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"


logger = logging.getLogger(__name__)


def get_binding_energies(
    mol_operation_type,
    mol,
    index,
    bond_order=1,
    db=None,
    name="binding_energy_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    **kwargs
):
    # TODO: test with different charges and spin multiplicities when
    #  deriving molecules
    # TODO: include an option to use free energy instead of SCF energy
    # mol_operation_type = [], mol = [], index = [], skips = [[], []],
    # mol_name = []
    # order of the indices should be consistent with the order of the mols
    # process_mol_func applies to both input molecules, so if set to False,
    # should give mol_name to each molecule, and if set to True, both will
    # take MolFormula even if mol_name is given to either or both
    fws = []
    labels = []
    working_dir = working_dir or os.getcwd()
    gout_keys = ["mol_1", "mol_2", "mol_linked"]
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    gaussian_inputs = handle_gaussian_inputs(
        {"opt": opt_gaussian_inputs, "freq": freq_gaussian_inputs},
        solvent_gaussian_inputs,
        solvent_properties,
    )
    opt_gaussian_inputs = gaussian_inputs["opt"]
    freq_gaussian_inputs = gaussian_inputs["freq"]

    if skips is None:
        skips = [None, None]
    check_result = []
    for i, j in enumerate(skips):
        if j:
            check_result.append(["final_energy"])
        else:
            check_result.append(None)

    parents = []
    for position, [operation, molecule, key, skip, check, molecule_name] in enumerate(
        zip(
            mol_operation_type,
            mol,
            gout_keys[:2],
            skips,
            check_result,
            kwargs.pop("mol_name", [None, None]),
        )
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
            check_result=check,
            **kwargs
        )
        fws += opt_freq_init_fws
        parents.append(len(fws))
        labels.append(label)

    final_mol_label = "{}_{}".format(labels[0], labels[1])

    kwargs.pop("process_mol_func", False)

    _, _, opt_freq_final_fws = common_fw(
        mol_operation_type="link_molecules",
        mol={
            "operation_type": ["get_from_run_dict", "get_from_run_dict"],
            "mol": gout_keys[:2],
            "index": index,
            "bond_order": bond_order,
        },
        working_dir=working_dir,
        db=db,
        filename=final_mol_label,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        process_mol_func=False,
        mol_name=final_mol_label,
        gout_key=gout_keys[-1],
        from_fw_spec=True,
        **kwargs
    )
    fws += opt_freq_final_fws
    links_dict = {fws[i - 1]: fws[-len(opt_freq_final_fws)] for i in parents}
    fw_analysis = Firework(
        BindingEnergytoDB(
            index=index,
            db=db,
            keys=gout_keys,
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            **{
                i: j
                for i, j in kwargs.items()
                if i
                in BindingEnergytoDB.required_params + BindingEnergytoDB.optional_params
            }
        ),
        parents=fws[:],
        name="{}-{}".format(final_mol_label, "binding_energy_analysis"),
        spec={"_launch_dir": os.path.join(working_dir, "analysis")},
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name="{}_{}".format(final_mol_label, name),
        links_dict=links_dict,
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS}
    )
