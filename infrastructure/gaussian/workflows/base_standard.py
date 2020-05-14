import os
import logging
import datetime

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
              db=None,
              process_mol_func=True,
              mol_name=None,
              dir_head=None,
              skip_opt_freq=False,
              **kwargs):
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
    # TODO: test with different charges and spin multiplicities when deriving molecules
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
                                          'analysis')})
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
              db=None,
              name="ip_ea_calculation",
              working_dir=None,
              opt_gaussian_inputs=None,
              freq_gaussian_inputs=None,
              cart_coords=True,
              skip_opt_freq=False,
              **kwargs):
    # TODO: check the working dirs ; they are wrong
    # TODO: check the links between the fireworks
    # TODO: check with the group if they want to skip opt and freq
    # TODO: check with the group if the order of optimizations is correct
    # TODO: cleanup the messy code
    # TODO: check with the group if one electron is exchanged only
    # TODO: check if the group if we want to skip either IP or EA

    fws = []
    parents = []
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    keys = ["neutral_gas", "anion_gas", "cation_gas",
            "neutral_sol", "anion_sol", "cation_sol"]

    if ref_charge != opt_gaussian_inputs.get("charge", ref_charge):
        raise Exception("The provided reference charge is not consistent with "
                        "the one found in the gaussian input parameters.")
    sol_opt_gaussian_inputs = opt_gaussian_inputs
    sol_freq_gaussian_inputs = freq_gaussian_inputs
    gaussian_inputs = [opt_gaussian_inputs, freq_gaussian_inputs,
                       sol_opt_gaussian_inputs, sol_freq_gaussian_inputs]
    for i in gaussian_inputs:
        i["charge"] = ref_charge
        i["spin_multiplicity"] = spin_multiplicities["neutral"]

    if "generic" in solvent_gaussian_inputs.lower() and not solvent_properties:
        raise Exception("A generic solvent is provided as an input without "
                        "specifying its parameters.")
    sol_opt_gaussian_inputs["route_parameters"]["SCRF"] = \
        solvent_gaussian_inputs
    sol_freq_gaussian_inputs["route_parameters"]["SCRF"] = \
        solvent_gaussian_inputs
    if solvent_properties:
        if "input_parameters" not in sol_opt_gaussian_inputs:
            sol_opt_gaussian_inputs["input_parameters"] = {}
        sol_opt_gaussian_inputs["input_parameters"].update(solvent_properties)

    mol, neutral_gas_fws = common_fw(mol_operation_type=mol_operation_type,
                                     mol=mol,
                                     working_dir=working_dir,
                                     db=db,
                                     opt_gaussian_inputs=opt_gaussian_inputs,
                                     freq_gaussian_inputs=freq_gaussian_inputs,
                                     cart_coords=cart_coords,
                                     oxidation_states=None,
                                     gout_key="neutral_gas",
                                     skip_opt_freq=skip_opt_freq,
                                     dir_structure=["Neutral", "Gas"],
                                     **kwargs)
    fws += neutral_gas_fws
    mol_formula = get_mol_formula(mol)
    # parents.append(len(fws))

    _, neutral_sol_fws = common_fw(mol_operation_type="get_from_run_dict",
                                   mol="neutral_gas",
                                   working_dir=working_dir,
                                   db=db,
                                   opt_gaussian_inputs=sol_opt_gaussian_inputs,
                                   freq_gaussian_inputs=sol_freq_gaussian_inputs,
                                   cart_coords=cart_coords,
                                   oxidation_states=None,
                                   process_mol_func=False,
                                   mol_name="{}_sol".format(mol_formula),
                                   dir_head=mol_formula,
                                   gout_key="neutral_sol",
                                   skip_opt_freq=False,
                                   dir_structure=["Neutral", "Solution"],
                                   from_fw_spec=True,
                                   **kwargs)
    fws += neutral_sol_fws
    # parents.append(len(fws))

    for i in gaussian_inputs[:2]:
        i["charge"] = ref_charge - num_electrons
        i["spin_multiplicity"] = spin_multiplicities["anion"]

    _, anion_gas_fws = common_fw(mol_operation_type="get_from_run_dict",
                                 mol="neutral_gas",
                                 working_dir=working_dir,
                                 db=db,
                                 opt_gaussian_inputs=opt_gaussian_inputs,
                                 freq_gaussian_inputs=freq_gaussian_inputs,
                                 cart_coords=cart_coords,
                                 oxidation_states=None,
                                 process_mol_func=False,
                                 mol_name="{}-_gas".format(mol_formula),
                                 dir_head=mol_formula,
                                 gout_key="anion_gas",
                                 skip_opt_freq=False,
                                 dir_structure=["Anion", "Gas"],
                                 from_fw_spec=True,
                                 **kwargs)
    fws += anion_gas_fws
    # parents.append(len(fws))

    for i in gaussian_inputs[:2]:
        i["charge"] = ref_charge + num_electrons
        i["spin_multiplicity"] = spin_multiplicities["cation"]

    _, cation_gas_fws = common_fw(mol_operation_type="get_from_run_dict",
                                  mol="neutral_gas",
                                  working_dir=working_dir,
                                  db=db,
                                  opt_gaussian_inputs=opt_gaussian_inputs,
                                  freq_gaussian_inputs=freq_gaussian_inputs,
                                  cart_coords=cart_coords,
                                  oxidation_states=None,
                                  process_mol_func=False,
                                  mol_name="{}+_gas".format(mol_formula),
                                  dir_head=mol_formula,
                                  gout_key="cation_gas",
                                  skip_opt_freq=False,
                                  dir_structure=["Cation", "Gas"],
                                  from_fw_spec=True,
                                  **kwargs)
    fws += cation_gas_fws
    # parents.append(len(fws))

    for i in gaussian_inputs[2:]:
        i["charge"] = ref_charge - num_electrons
        i["spin_multiplicity"] = spin_multiplicities["anion"]

    _, anion_sol_fws = common_fw(mol_operation_type="get_from_run_dict",
                                 mol="anion_gas",
                                 working_dir=working_dir,
                                 db=db,
                                 opt_gaussian_inputs=opt_gaussian_inputs,
                                 freq_gaussian_inputs=freq_gaussian_inputs,
                                 cart_coords=cart_coords,
                                 oxidation_states=None,
                                 process_mol_func=False,
                                 mol_name="{}-_sol".format(mol_formula),
                                 dir_head=mol_formula,
                                 gout_key="anion_sol",
                                 skip_opt_freq=False,
                                 dir_structure=["Anion", "Solution"],
                                 from_fw_spec=True,
                                 **kwargs)
    fws += anion_sol_fws
    # parents.append(len(fws))

    for i in gaussian_inputs[2:]:
        i["charge"] = ref_charge + num_electrons
        i["spin_multiplicity"] = spin_multiplicities["cation"]

    _, cation_sol_fws = common_fw(mol_operation_type="get_from_run_dict",
                                  mol="cation_gas",
                                  working_dir=working_dir,
                                  db=db,
                                  opt_gaussian_inputs=opt_gaussian_inputs,
                                  freq_gaussian_inputs=freq_gaussian_inputs,
                                  cart_coords=cart_coords,
                                  oxidation_states=None,
                                  process_mol_func=False,
                                  mol_name="{}+_sol".format(mol_formula),
                                  dir_head=mol_formula,
                                  gout_key="cation_sol",
                                  skip_opt_freq=False,
                                  dir_structure=["Cation", "Solution"],
                                  from_fw_spec=True,
                                  **kwargs)
    fws += cation_sol_fws
    # parents.append(len(fws))
    links_dict = {fws[1]: [fws[2], fws[4], fws[6]], fws[5]: fws[8],
                  fws[7]: fws[10]}
    fw_analysis = Firework(
        IPEAtoDB(num_electrons=num_electrons,
                 db=db,
                 **{i: j for i, j in kwargs.items()
                    if i in BindingEnergytoDB.required_params +
                    BindingEnergytoDB.optional_params}),
        parents=fws[:],
        name="{}-{}".format(mol_formula,
                            "ip_ea_analysis"),
        spec={'_launch_dir': os.path.join(working_dir, mol_formula, "analysis")})
    fws.append(fw_analysis)

    return Workflow(fws,
                    name="{}_{}".format(mol_formula, name),
                    links_dict=links_dict,
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})
