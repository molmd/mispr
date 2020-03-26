import os
import logging

from fireworks import Firework, Workflow
from pymatgen.core.structure import Molecule
from infrastructure.gaussian.firetasks.geo_transformation import \
    BindingEnergytoDB
from infrastructure.gaussian.firetasks.parse_outputs import ProcessRun
from infrastructure.gaussian.fireworks.core_standard import CalcFromMolFW, \
    CalcFromRunsDBFW
from infrastructure.gaussian.utils.utils import process_mol, get_job_name, \
    get_mol_formula, process_run

logger = logging.getLogger(__name__)

WORKFLOW_KWARGS = ["links_dict", "name", "metadata", "created_on",
                   "updated_on", "fw_states"]

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
                               working_dir=working_dir,
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
                                   working_dir=working_dir,
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
        # mol_operation_type = get_from_file, get_from_gout, get_from_run_id,
        # get_from_run query, get_from_gout, get_from_run_dict only
        # I'm repeating the same thing twice here!
        # TODO: raise error if the mol_operation_type is not from above
        run = process_run(mol_operation_type, mol, db=db,
                          working_dir=working_dir)
        mol = process_mol("get_from_run_dict", run, db=db,
                          working_dir=working_dir)
        spec = kwargs.pop('spec', {})
        if "tag" in kwargs:
            spec.update({'tag': kwargs["tag"]})
        print(kwargs)
        fws.append(Firework(ProcessRun(run=run,
                                       operation_type="get_from_run_dict",
                                       db=db,
                                       working_dir=working_dir,
                                       **{i: j for i, j in kwargs.items() if
                                          i in
                                          ProcessRun.required_params +
                                          ProcessRun.optional_params}),
                            name=get_job_name(mol, "process_run"),
                            spec=spec)
                   )
        print(kwargs)
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
                    save_to_db=False,
                    update_duplicates=False,
                    oxidation_states=None,
                    skip_opt_freq=False,
                    **kwargs):
    print(kwargs)
    fws = []
    working_dir = working_dir or os.getcwd()

    mol, opt_freq_fws = common_fw(mol_operation_type=mol_operation_type,
                                  mol=mol,
                                  working_dir=working_dir,
                                  db=db,
                                  opt_gaussian_inputs=opt_gaussian_inputs,
                                  freq_gaussian_inputs=freq_gaussian_inputs,
                                  cart_coords=cart_coords,
                                  save_to_db=save_to_db,
                                  update_duplicates=update_duplicates,
                                  oxidation_states=oxidation_states,
                                  skip_opt_freq=skip_opt_freq,
                                  **kwargs)
    print(kwargs)
    fws += opt_freq_fws
    esp_gaussian_inputs = esp_gaussian_inputs or {}
    if "route_parameters" not in esp_gaussian_inputs:
        esp_gaussian_inputs.update({"route_parameters": {"pop": "MK",
                                                         "iop(6/50=1)": None}})
    # input_parameters from a previous run are overwritten
    if "input_parameters" not in esp_gaussian_inputs:
        esp_gaussian_inputs.update({"input_parameters": {"molesp": None}})
    mol_formula = get_mol_formula(mol)
    print(kwargs)
    esp_fw = CalcFromRunsDBFW(db,
                              input_file="{}_esp.com".format(mol_formula),
                              output_file="{}_esp.out".format(mol_formula),
                              name=get_job_name(mol, "esp"),
                              parents=fws[:],
                              gaussian_input_params=esp_gaussian_inputs,
                              working_dir=working_dir,
                              cart_coords=cart_coords,
                              spec={"proceed": {"has_gaussian_completed": True,
                                                "stationary_type": "Minimum"}},
                              **kwargs
                              )
    fws.append(esp_fw)
    return Workflow(fws,
                    name=get_job_name(mol, name),
                    **kwargs)


def get_nmr_tensors(mol_operation_type,
                    mol,
                    db=None,
                    name="nmr_tensor_calculation",
                    working_dir=None,
                    opt_gaussian_inputs=None,
                    freq_gaussian_inputs=None,
                    nmr_gaussian_inputs=None,
                    cart_coords=True,
                    save_to_db=False,
                    update_duplicates=False,
                    oxidation_states=None,
                    skip_opt_freq=False,
                    **kwargs):
    fws = []
    working_dir = working_dir or os.getcwd()

    mol, opt_freq_fws = common_fw(mol_operation_type=mol_operation_type,
                                  mol=mol,
                                  working_dir=working_dir,
                                  db=db,
                                  opt_gaussian_inputs=opt_gaussian_inputs,
                                  freq_gaussian_inputs=freq_gaussian_inputs,
                                  cart_coords=cart_coords,
                                  save_to_db=save_to_db,
                                  update_duplicates=update_duplicates,
                                  oxidation_states=oxidation_states,
                                  skip_opt_freq=skip_opt_freq,
                                  **kwargs)
    fws += opt_freq_fws
    nmr_gaussian_inputs = nmr_gaussian_inputs or {}
    if "route_parameters" not in nmr_gaussian_inputs:
        nmr_gaussian_inputs.update({"route_parameters": {"NMR": "GIAO"}})
    mol_formula = get_mol_formula(mol)

    nmr_fw = CalcFromRunsDBFW(db,
                              input_file="{}_nmr.com".format(mol_formula),
                              output_file="{}_nmr.out".format(mol_formula),
                              name=get_job_name(mol, "nmr"),
                              parents=fws[1],
                              gaussian_input_params=nmr_gaussian_inputs,
                              working_dir=working_dir,
                              cart_coords=cart_coords,
                              spec={"proceed": {"has_gaussian_completed": True,
                                                "stationary_type": "Minimum"}},
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
                         save_to_db=False,
                         update_duplicates=False,
                         oxidation_states=None,
                         skip_opt_freq=None,
                         **kwargs):
    # mol_operation_type = [], mol = [], index = []
    # order of the indices should be consistent with the order of the mols

    if skip_opt_freq is None:
        skip_opt_freq = [False, False]
    fws = []
    molecules = []
    working_dir = working_dir or os.getcwd()
    key1 = 'linking_mol'
    key2 = 'linked_mol'
    key3 = 'final_mol'

    for position, [operation, molecule, key, skip] in \
            enumerate(
                zip(mol_operation_type, mol, [key1, key2], skip_opt_freq)):
        mol_object, list_fws_1 = common_fw(mol_operation_type=operation,
                                           mol=molecule,
                                           working_dir=working_dir,
                                           db=db,
                                           opt_gaussian_inputs=opt_gaussian_inputs,
                                           freq_gaussian_inputs=freq_gaussian_inputs,
                                           cart_coords=cart_coords,
                                           save_to_db=save_to_db,
                                           update_duplicates=update_duplicates,
                                           oxidation_states=oxidation_states,
                                           gout_key=key,
                                           skip_opt_freq=skip,
                                           **kwargs)

        fws += list_fws_1
        molecules.append(mol_object)

    final_mol_name = "{}_{}".format(get_mol_formula(molecules[0]),
                                    get_mol_formula(molecules[1]))

    _, list_fws_2 = common_fw(mol_operation_type="link_molecules",
                              mol={"operation_type": ["get_from_run_dict",
                                                      "get_from_run_dict"],
                                   "mol": [key1, key2],
                                   "index": index,
                                   "bond_order": bond_order},
                              working_dir=working_dir,
                              db=db,
                              opt_gaussian_inputs=opt_gaussian_inputs,
                              freq_gaussian_inputs=freq_gaussian_inputs,
                              cart_coords=cart_coords,
                              save_to_db=save_to_db,
                              update_duplicates=update_duplicates,
                              oxidation_states=oxidation_states,
                              process_mol_func=False,
                              mol_name=final_mol_name,
                              gout_key=key3,
                              from_fw_spec=True)

    fws += list_fws_2

    fw_analysis = Firework(BindingEnergytoDB(keys=[key1, key2, key3],
                                             prop="final_energy",
                                             main_run_key=key3,
                                             new_prop="binding_energy_{}_{}_eV".
                                             format(
                                                 molecules[0].species[index[0]],
                                                 molecules[1].species[
                                                     index[1]]),
                                             db=db),
                           parents=fws[:],
                           name="{}-{}".format(final_mol_name,
                                               "binding_energy_analysis"))
    fws.append(fw_analysis)

    name = "{}_{}".format(final_mol_name, name)
    links_dict = {fws[1]: fws[4], fws[3]: fws[4]}

    return Workflow(fws,
                    name=name,
                    links_dict=links_dict,
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})
