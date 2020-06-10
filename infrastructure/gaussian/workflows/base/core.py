import os
import logging
import itertools
from copy import deepcopy

from fireworks import Firework, Workflow
from infrastructure.gaussian.firetasks.parse_outputs import \
    BindingEnergytoDB, IPEAtoDB
from infrastructure.gaussian.firetasks.parse_outputs import ProcessRun
from infrastructure.gaussian.fireworks.core_standard import CalcFromMolFW, \
    CalcFromRunsDBFW
from infrastructure.gaussian.utils.utils import process_mol, get_job_name, \
    get_mol_formula, process_run, recursive_relative_to_absolute_path

logger = logging.getLogger(__name__)

WORKFLOW_KWARGS = Workflow.__init__.__code__.co_varnames

STANDARD_OPT_GUASSIAN_INPUT = {"functional": "B3LYP",
                               "basis_set": "6-31G(d)",
                               "route_parameters": {"Opt": None},
                               "link0_parameters": {"%chk": "checkpoint.chk",
                                                    "%mem": "1000MW",
                                                    "%NProcShared": "24"}}


def common_fw(mol_operation_type,
              mol,
              working_dir,
              opt_gaussian_inputs,
              freq_gaussian_inputs,
              cart_coords,
              oxidation_states,
              gout_key=None,
              db=None,
              process_mol_func=True,
              mol_name=None,
              dir_head=None,
              skip_opt_freq=False,
              **kwargs):
    # TODO: checking criteria not working
    fws = []
    if not skip_opt_freq:
        if process_mol_func:
            mol = process_mol(mol_operation_type, mol, db=db,
                              working_dir=working_dir)
            mol_operation_type = 'get_from_mol'
            mol_formula = get_mol_formula(mol)
            opt_job_name = get_job_name(mol, "optimization")
            freq_job_name = get_job_name(mol, "frequency")
            label = mol_formula
        elif mol_name:
            opt_job_name = "{}_optimization".format(mol_name)
            freq_job_name = "{}_frequency".format(mol_name)
            label = mol_name
        else:
            opt_job_name = "optimization"
            freq_job_name = "frequency"
            label = "mol"
        if not dir_head:
            dir_head = label
        dir_struct = [dir_head] + kwargs.get('dir_structure', [])
        working_dir = os.path.join(working_dir, *dir_struct)
        opt_input_file = f"{label}_opt.com"
        opt_output_file = f"{label}_opt.out"
        freq_input_file = f"{label}_freq.com"
        freq_output_file = f"{label}_freq.out"

        opt_gaussian_inputs = opt_gaussian_inputs or {}
        opt_gaussian_inputs = {**STANDARD_OPT_GUASSIAN_INPUT,
                               **opt_gaussian_inputs}
        if "opt" not in [i.lower() for i in
                         opt_gaussian_inputs["route_parameters"]]:
            raise ValueError("The Opt keyword is missing from the input file")
        opt_fw = CalcFromMolFW(mol=mol,
                               mol_operation_type=mol_operation_type,
                               db=db,
                               name=opt_job_name,
                               working_dir=os.path.join(working_dir,
                                                        "Optimization"),
                               input_file=opt_input_file,
                               output_file=opt_output_file,
                               gaussian_input_params=opt_gaussian_inputs,
                               cart_coords=cart_coords,
                               oxidation_states=oxidation_states,
                               gout_key=gout_key+"_opt",
                               **kwargs
                               )
        fws.append(opt_fw)

        # if no freq_gaussian_inputs are provided, parameters from prev opt
        # are used except for the route parameters which are replaced with
        # the Freq keyword
        freq_gaussian_inputs = freq_gaussian_inputs or {}
        if "route_parameters" not in freq_gaussian_inputs:
            freq_gaussian_inputs.update({"route_parameters": {"Freq": None}})
        if "freq" not in [i.lower() for i in
                          freq_gaussian_inputs["route_parameters"]]:
            raise ValueError('The Freq keyword is missing from the input file')

        freq_fw = CalcFromRunsDBFW(db=db,
                                   name=freq_job_name,
                                   parents=opt_fw,
                                   gaussian_input_params=freq_gaussian_inputs,
                                   working_dir=os.path.join(working_dir,
                                                            "Frequency"),
                                   input_file=freq_input_file,
                                   output_file=freq_output_file,
                                   cart_coords=cart_coords,
                                   gout_key=gout_key,
                                   spec={
                                       "proceed": {
                                           "has_gaussian_completed": True}},
                                   **kwargs
                                   )
        fws.append(freq_fw)

    else:
        if mol_operation_type not in ["get_from_gout", "get_from_file",
                                      "get_from_run_dict", "get_from_run_id",
                                      "get_from_run_query"]:
            raise ValueError("no Gaussian output provided; to skip "
                             "optimization and frequency, you need to input "
                             "the molecule in any of the supported Gaussian "
                             "output formats")
        else:
            run = process_run(mol_operation_type, mol, db=db,
                              working_dir=working_dir)
            mol = process_mol("get_from_run_dict", run, db=db,
                              working_dir=working_dir)
            spec = kwargs.pop('spec', {})
            label = get_mol_formula(mol)
            if not dir_head:
                dir_head = label
            dir_struct = [dir_head] + kwargs.get('dir_structure', [])
            working_dir = os.path.join(working_dir, *dir_struct)
            if "tag" in kwargs:
                spec.update({'tag': kwargs["tag"]})
            spec.update({'_launch_dir': working_dir})
            fws.append(Firework(ProcessRun(run=run,
                                           operation_type="get_from_run_dict",
                                           db=db,
                                           gout_key=gout_key,
                                           **{i: j for i, j in kwargs.items() if
                                              i in
                                              ProcessRun.required_params +
                                              ProcessRun.optional_params}),
                                name=get_job_name(mol, "process_run"),
                                spec=spec)
                       )
    return mol, fws


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

    mol, opt_freq_fws = common_fw(mol_operation_type=mol_operation_type,
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
    mol_formula = get_mol_formula(mol)
    esp_gaussian_inputs = esp_gaussian_inputs or {}
    if "route_parameters" not in esp_gaussian_inputs:
        esp_gaussian_inputs.update({"route_parameters": {"pop": "MK",
                                                         "iop(6/50=1)": None}})
    # input_parameters from a previous run are overwritten
    if "input_parameters" not in esp_gaussian_inputs:
        mol_esp = os.path.join(
            working_dir, "{}_esp".format(
                os.path.join(working_dir, mol_formula, "ESP", mol_formula)))
        esp_gaussian_inputs.update({"input_parameters": {mol_esp: None}})

    if not skip_opt_freq:
        spec = {"proceed": {"has_gaussian_completed": True,
                            "stationary_type": "Minimum"}}
    else:
        spec = {"proceed": {"has_gaussian_completed": True}}

    esp_fw = CalcFromRunsDBFW(db,
                              input_file="{}_esp.com".format(mol_formula),
                              output_file="{}_esp.out".format(mol_formula),
                              name=get_job_name(mol, "esp"),
                              parents=fws[:],
                              gaussian_input_params=esp_gaussian_inputs,
                              working_dir=os.path.join(working_dir, mol_formula,
                                                       "ESP"),
                              cart_coords=cart_coords,
                              spec=spec,
                              **kwargs
                              )
    fws.append(esp_fw)
    return Workflow(fws,
                    name=get_job_name(mol, name),
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})


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

    mol, opt_freq_fws = common_fw(mol_operation_type=mol_operation_type,
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
    mol_formula = get_mol_formula(mol)
    nmr_gaussian_inputs = nmr_gaussian_inputs or {}
    if "route_parameters" not in nmr_gaussian_inputs:
        nmr_gaussian_inputs.update({"route_parameters": {"NMR": "GIAO"}})

    if not skip_opt_freq:
        spec = {"proceed": {"has_gaussian_completed": True,
                            "stationary_type": "Minimum"}}
    else:
        spec = {"proceed": {"has_gaussian_completed": True}}

    nmr_fw = CalcFromRunsDBFW(db,
                              input_file="{}_nmr.com".format(mol_formula),
                              output_file="{}_nmr.out".format(mol_formula),
                              name=get_job_name(mol, "nmr"),
                              parents=fws[:],
                              gaussian_input_params=nmr_gaussian_inputs,
                              working_dir=os.path.join(working_dir, mol_formula,
                                                       "NMR"),
                              cart_coords=cart_coords,
                              spec=spec,
                              **kwargs
                              )
    fws.append(nmr_fw)
    return Workflow(fws,
                    name=get_job_name(mol, name),
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})


def get_binding_energies(mol_operation_type,
                         mol,
                         index,
                         bond_order=1,
                         db=None,
                         name="binding_energy_calculation",
                         working_dir=None,
                         opt_gaussian_inputs=None,
                         freq_gaussian_inputs=None,
                         cart_coords=True,
                         oxidation_states=None,
                         skip_opt_freq=None,
                         **kwargs):
    # TODO: test with different charges and spin multiplicities when
    #  deriving molecules
    # mol_operation_type = [], mol = [], index = []
    # order of the indices should be consistent with the order of the mols
    fws = []
    molecules = []
    working_dir = working_dir or os.getcwd()
    keys = ["mol_1", "mol_2", "mol_linked"]
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    if skip_opt_freq is None:
        skip_opt_freq = [False, False]
    parents = []
    for position, [operation, molecule, key, skip] in \
            enumerate(
                zip(mol_operation_type, mol, keys[:2], skip_opt_freq)):
        mol_object, opt_freq_init_fws = \
            common_fw(mol_operation_type=operation,
                      mol=molecule,
                      working_dir=working_dir,
                      db=db,
                      opt_gaussian_inputs=opt_gaussian_inputs,
                      freq_gaussian_inputs=freq_gaussian_inputs,
                      cart_coords=cart_coords,
                      oxidation_states=oxidation_states,
                      gout_key=key,
                      skip_opt_freq=skip,
                      **kwargs)
        fws += opt_freq_init_fws
        parents.append(len(fws))
        molecules.append(mol_object)

    final_mol_formula = "{}_{}".format(get_mol_formula(molecules[0]),
                                       get_mol_formula(molecules[1]))

    _, opt_freq_final_fws = common_fw(mol_operation_type="link_molecules",
                                      mol={"operation_type": [
                                          "get_from_run_dict",
                                          "get_from_run_dict"],
                                          "mol": keys[:2],
                                          "index": index,
                                          "bond_order": bond_order},
                                      working_dir=working_dir,
                                      db=db,
                                      filename=final_mol_formula,
                                      opt_gaussian_inputs=opt_gaussian_inputs,
                                      freq_gaussian_inputs=freq_gaussian_inputs,
                                      cart_coords=cart_coords,
                                      oxidation_states=oxidation_states,
                                      process_mol_func=False,
                                      mol_name=final_mol_formula,
                                      gout_key=keys[-1],
                                      from_fw_spec=True,
                                      **kwargs)
    fws += opt_freq_final_fws
    links_dict = {fws[i - 1]: fws[-len(opt_freq_final_fws)] for i in parents}
    fw_analysis = Firework(
        BindingEnergytoDB(index=index,
                          db=db,
                          **{i: j for i, j in kwargs.items()
                             if i in BindingEnergytoDB.required_params +
                             BindingEnergytoDB.optional_params}),
        parents=fws[:],
        name="{}-{}".format(final_mol_formula,
                            "binding_energy_analysis"),
        spec={'_launch_dir': os.path.join(working_dir, final_mol_formula,
                                          'Analysis')})
    fws.append(fw_analysis)

    return Workflow(fws,
                    name="{}_{}".format(final_mol_formula, name),
                    links_dict=links_dict,
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})


def get_ip_ea(mol_operation_type,
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
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    if ref_charge != opt_gaussian_inputs.get("charge", ref_charge):
        raise Exception("The provided reference charge is not consistent with "
                        "the one found in the gaussian input parameters.")

    # TODO: Cleanup this part
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
        fws += _create_fireworks(state, phase)
    mol_formula = get_mol_formula(mol)
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
                    if i in BindingEnergytoDB.required_params +
                    BindingEnergytoDB.optional_params}),
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
