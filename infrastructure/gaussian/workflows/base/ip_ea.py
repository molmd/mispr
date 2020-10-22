import os
import logging
import itertools

from copy import deepcopy

from fireworks import Firework, Workflow

from infrastructure.gaussian.utils.utils import get_mol_formula, \
    recursive_relative_to_absolute_path
from infrastructure.gaussian.firetasks.parse_outputs import IPEAtoDB
from infrastructure.gaussian.workflows.base.core import common_fw, \
    WORKFLOW_KWARGS

logger = logging.getLogger(__name__)


def _recursive_create_fireworks(num_steps):
    if num_steps == 1:
        pass


def get_ip_ea_single_step(mol_operation_type,
                          mol,
                          ref_charge,
                          spin_multiplicities,
                          num_electrons=1,
                          solvent_gaussian_inputs=None,
                          solvent_properties=None,
                          states=None,
                          phases=None,
                          electrode_potentials=None,
                          db=None,
                          name="ip_ea_calculation",
                          working_dir=None,
                          opt_gaussian_inputs=None,
                          freq_gaussian_inputs=None,
                          cart_coords=True,
                          skip_opt_freq=False,
                          **kwargs):
    def _create_fireworks(local_state, local_phase):
        local_key = "{}_{}".format(local_state.lower(), local_phase.lower())
        parent_state, parent_phase = _get_parent(local_state, local_phase)
        parent_key = None
        if parent_state in states and parent_phase in phases:
            parent_key = "{}_{}".format(parent_state.lower(),
                                        parent_phase.lower())
        if not parent_key:
            local_mol_op_type = mol_operation_type
            local_mol = mol
            local_skip = skip_opt_freq
            local_process_mol_func = True
            local_mol_name = None
            local_from_fw_spec = False
            local_dir_head = None
        else:
            local_mol_op_type = 'get_from_run_dict'
            local_mol = parent_key
            local_skip = False
            local_process_mol_func = False
            local_mol_name = "{}_{}".format(mol_formula, local_phase.lower())
            local_from_fw_spec = True
            local_dir_head = mol_formula

        add_charge = 1 if local_state.lower() == 'cation' \
            else -1 if local_state.lower() == 'anion' else 0
        add_charge *= num_electrons

        local_opt_gins = deepcopy(opt_gaussian_inputs)
        local_freq_gins = deepcopy(freq_gaussian_inputs)
        local_opt_gins['charge'] = ref_charge + add_charge
        local_freq_gins['charge'] = ref_charge + add_charge
        local_opt_gins["spin_multiplicity"] = spin_multiplicities[local_state]
        local_freq_gins["spin_multiplicity"] = spin_multiplicities[local_state]

        if local_phase.lower() == 'solution':
            if "generic" in solvent_gaussian_inputs.lower() \
                    and not solvent_properties:
                raise Exception(
                    "A generic solvent is provided as an input without "
                    "specifying its parameters.")
            local_opt_gins["route_parameters"]["SCRF"] = \
                solvent_gaussian_inputs
            local_freq_gins["route_parameters"]["SCRF"] = \
                solvent_gaussian_inputs
            if solvent_properties:
                if "input_parameters" not in local_opt_gins:
                    local_opt_gins["input_parameters"] = {}
                local_opt_gins["input_parameters"].update(solvent_properties)

        local_mol, local_fws = \
            common_fw(mol_operation_type=local_mol_op_type,
                      mol=local_mol,
                      working_dir=working_dir,
                      db=db,
                      opt_gaussian_inputs=local_opt_gins,
                      freq_gaussian_inputs=local_freq_gins,
                      cart_coords=cart_coords,
                      oxidation_states=None,
                      process_mol_func=local_process_mol_func,
                      mol_name=local_mol_name,
                      dir_head=local_dir_head,
                      gout_key=local_key,
                      skip_opt_freq=local_skip,
                      dir_structure=[local_state, local_phase],
                      from_fw_spec=local_from_fw_spec,
                      **kwargs)
        fireworks_dict[local_key] = local_fws
        if parent_key:
            parents_dict[local_key] = parent_key
            return local_fws  # local_mol, local_fws, fireworks_dict
        else:
            return local_mol, local_fws

    def _get_parent(local_state, local_phase):
        if local_phase.lower() == "solution":
            if 'gas' in [k.lower() for k in states]:
                return local_state, "gas"
        if local_state.lower() == "cation" or local_state.lower() == "anion":
            return "reference", local_state
        else:
            return None, None

    fws = []
    fireworks_dict = {}
    links_dict = {}
    parents_dict = {}
    mol_object = None
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    if ref_charge != opt_gaussian_inputs.get("charge", ref_charge):
        raise Exception("The provided reference charge is not consistent with "
                        "the one found in the gaussian input parameters.")

    if states is None:
        states = ["reference", "cation", "anion"]
    if phases is None:
        phases = ["gas", "solution"]

    for state in states:
        if state.lower() not in ["reference", "cation", "anion"]:
            raise ValueError("The provided states are not supported. Supported "
                             "ones are reference, cation, and/or anion.")
    for phase in phases:
        if phase.lower() not in ["gas", "solution"]:
            raise ValueError("The provided phases are not supported. Supported "
                             "ones are gas and/or solution.")

    for state, phase in itertools.product(states, phases):
        cr_out = _create_fireworks(state, phase)
        if len(cr_out) == 2:
            mol_object, opt_freq_fws = cr_out
        else:
            opt_freq_fws = cr_out
        fws += opt_freq_fws
    mol_formula = get_mol_formula(mol_object)
    for i, j in parents_dict.items():
        links_dict[fireworks_dict[j][1]] = \
            links_dict.get(fireworks_dict[j][1], []) + [fireworks_dict[j][0]]
    fw_analysis = Firework(
        IPEAtoDB(num_electrons=num_electrons,
                 states=states,
                 phases=phases,
                 solvent_gaussian_inputs=solvent_gaussian_inputs,
                 solvent_properties=solvent_properties,
                 electrode_potentials=electrode_potentials,
                 parents_dict=parents_dict,
                 db=db,
                 **{i: j for i, j in kwargs.items()
                    if i in IPEAtoDB.required_params +
                    IPEAtoDB.optional_params}),
        parents=fws[:],
        name="{}-{}".format(mol_formula,
                            "ip_ea_analysis"),
        spec={
            '_launch_dir': os.path.join(working_dir, mol_formula, "Analysis")})
    fws.append(fw_analysis)

    return Workflow(fws,
                    name="{}_{}".format(mol_formula, name),
                    links_dict=links_dict,
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})


def get_ip_ea_multi_step(mol_operation_type,
                         mol,
                         ref_charge,
                         spin_multiplicities,
                         num_electrons=1,
                         solvent_gaussian_inputs=None,
                         solvent_properties=None,
                         states=None,
                         phases=None,
                         electrode_potentials=None,
                         db=None,
                         name="ip_ea_calculation",
                         working_dir=None,
                         opt_gaussian_inputs=None,
                         freq_gaussian_inputs=None,
                         cart_coords=True,
                         skip_opt_freq=False,
                         **kwargs):
    for i in range(num_electrons):
        local_num_electrons = i + 1
        wf = get_ip_ea_single_step(mol_operation_type=mol_operation_type,
                                   mol=mol,
                                   ref_charge=ref_charge,
                                   spin_multiplicities=spin_multiplicities,
                                   num_electrons=1,
                                   solvent_gaussian_inputs=solvent_gaussian_inputs,
                                   solvent_properties=solvent_properties,
                                   states=states,
                                   phases=phases,
                                   electrode_potentials=electrode_potentials,
                                   db=db,
                                   name=name,
                                   working_dir=working_dir,
                                   opt_gaussian_inputs=opt_gaussian_inputs,
                                   freq_gaussian_inputs=freq_gaussian_inputs,
                                   cart_coords=cart_coords,
                                   skip_opt_freq=skip_opt_freq,
                                   **kwargs)
    pass
