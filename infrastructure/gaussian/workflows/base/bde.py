import os

from fireworks import Firework, Workflow

from infrastructure.gaussian.utils.utils import \
    recursive_relative_to_absolute_path, check_solvent_inputs, \
    add_solvent_inputs, get_job_name
from infrastructure.gaussian.firetasks.parse_outputs import BDEtoDB
from infrastructure.gaussian.fireworks.core_standard import BreakMolFW
from infrastructure.gaussian.workflows.base.core import common_fw, \
    WORKFLOW_KWARGS


def get_bde(mol_operation_type,
            mol,
            ref_charge=0,
            fragment_charges=None,
            bonds=None,
            db=None,
            name="bde_calculation",
            working_dir=None,
            opt_gaussian_inputs=None,
            freq_gaussian_inputs=None,
            solvent_gaussian_inputs=None,
            solvent_properties=None,
            cart_coords=True,
            oxidation_states=None,
            skip_opt_freq=False,
            **kwargs):
    # optional breaking bonds in subsequent species? NO
    fws = []
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)
    gout_key = "ref_mol"
    gins = [opt_gaussian_inputs, freq_gaussian_inputs]
    check_solvent_inputs(gins)

    if solvent_gaussian_inputs:
        opt_gaussian_inputs, freq_gaussian_inputs = \
            add_solvent_inputs(gins,
                               solvent_gaussian_inputs,
                               solvent_properties)

    _, label, opt_freq_fws = common_fw(mol_operation_type=mol_operation_type,
                                       mol=mol,
                                       charge=ref_charge,
                                       working_dir=working_dir,
                                       dir_structure=["principle_mol"],
                                       db=db,
                                       opt_gaussian_inputs=opt_gaussian_inputs,
                                       freq_gaussian_inputs=freq_gaussian_inputs,
                                       cart_coords=cart_coords,
                                       oxidation_states=oxidation_states,
                                       skip_opt_freq=skip_opt_freq,
                                       gout_key=gout_key,
                                       **kwargs)
    fws += opt_freq_fws

    break_fw = BreakMolFW(mol=gout_key,
                          mol_operation_type="get_from_run_dict",
                          from_fw_spec=True,
                          bonds=bonds,
                          ref_charge=ref_charge,
                          fragment_charges=fragment_charges,
                          db=db,
                          calc_frags=True,
                          opt_gaussian_inputs=opt_gaussian_inputs,
                          freq_gaussian_inputs=freq_gaussian_inputs,
                          cart_coords=cart_coords,
                          name=get_job_name(label, "breaking"),
                          parents=fws[:],
                          working_dir=os.path.join(working_dir,
                                                   label,
                                                   "fragments"),
                          **kwargs)
    fws.append(break_fw)

    fw_analysis = Firework(
        BDEtoDB(principle_mol_key=gout_key,
                db=db,
                solvent_gaussian_inputs=solvent_gaussian_inputs,
                solvent_properties=solvent_properties,
                **{i: j for i, j in kwargs.items()
                   if i in BDEtoDB.required_params +
                   BDEtoDB.optional_params}),
        parents=fws[:],
        name="{}-{}".format(label, "bde_analysis"),
        spec={'_launch_dir': os.path.join(working_dir, 'analysis')})
    fws.append(fw_analysis)

    return Workflow(fws,
                    name="{}_{}".format(label, name),
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})
