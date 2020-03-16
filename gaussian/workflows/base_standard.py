import os
import logging

from fireworks import Firework, Workflow
from infrastructure.gaussian.fireworks.core_standard import CalcFromMolFW, \
    CalcFromRunsDBFW
from infrastructure.gaussian.utils.utils import process_mol, get_job_name, \
    get_mol_formula

logger = logging.getLogger(__name__)

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
              **kwargs):
    if process_mol_func:
        mol = process_mol(mol_operation_type, mol, db=db,
                          working_dir=working_dir)
        mol_operation_type = 'get_from_mol'
        mol_formula = get_mol_formula(mol)
        opt_job_name = get_job_name(mol, "optimization")
        freq_job_name = get_job_name(mol, "frequency")
        opt_input_file = "{}_opt.com".format(mol_formula)
        opt_output_file = "{}_opt.out".format(mol_formula)
        freq_input_file = "{}_freq.com".format(mol_formula)
        freq_output_file = "{}_freq.out".format(mol_formula)
    else:
        opt_job_name = "optimization"
        freq_job_name = "frequency"
        opt_input_file = "mol_opt.com"
        opt_output_file = "mol_opt.out"
        freq_input_file = "mol_freq.com"
        freq_output_file = "mol_freq.out"
    print(mol_operation_type, mol)
    opt_gaussian_inputs = opt_gaussian_inputs or {}
    opt_gaussian_inputs = {**STANDARD_OPT_GUASSIAN_INPUT, **opt_gaussian_inputs}
    if "opt" not in [i.lower() for i in
                     opt_gaussian_inputs["route_parameters"]]:
        raise ValueError("The Opt keyword is missing from the input file")
    fw1 = CalcFromMolFW(mol=mol,
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

    # if no freq_gaussian_inputs are provided, parameters from prev opt are used
    # except for the route parameters which are replaced with the Freq keyword
    freq_gaussian_inputs = freq_gaussian_inputs or {}
    if "route_parameters" not in freq_gaussian_inputs:
        freq_gaussian_inputs.update({"route_parameters": {"Freq": None}})
    if "freq" not in [i.lower() for i in
                      freq_gaussian_inputs["route_parameters"]]:
        raise ValueError('The Freq keyword is missing from the input file')

    fw2 = CalcFromRunsDBFW(db=db,
                           name=freq_job_name,
                           parents=fw1,
                           gaussian_input_params=freq_gaussian_inputs,
                           working_dir=working_dir,
                           input_file=freq_input_file,
                           output_file=freq_output_file,
                           cart_coords=cart_coords,
                           spec={
                               "proceed": {"has_gaussian_completed": True}},
                           **kwargs
                           )
    return mol, [fw1, fw2]


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
                    save_mol_file=False,
                    oxidation_states=None,
                    **kwargs):
    fws = []
    working_dir = working_dir or os.getcwd()

    mol, list_fws = common_fw(mol_operation_type=mol_operation_type,
                              mol=mol,
                              working_dir=working_dir,
                              db=db,
                              opt_gaussian_inputs=opt_gaussian_inputs,
                              freq_gaussian_inputs=freq_gaussian_inputs,
                              cart_coords=cart_coords,
                              save_to_db=save_to_db,
                              update_duplicates=update_duplicates,
                              oxidation_states=oxidation_states,
                              save_mol_file=save_mol_file,
                              **kwargs)
    fws += list_fws
    esp_gaussian_inputs = esp_gaussian_inputs or {}
    if "route_parameters" not in esp_gaussian_inputs:
        esp_gaussian_inputs.update({"route_parameters": {"pop": "MK",
                                                         "iop(6/50=1)": None}})
    # input_parameters from a previous run are overwritten
    if "input_parameters" not in esp_gaussian_inputs:
        esp_gaussian_inputs.update({"input_parameters": {"molesp": None}})
    fws.append(
        CalcFromRunsDBFW(db,
                         input_file="mol_esp.com",
                         output_file="mol_esp.out",
                         name=get_job_name(mol, "esp"),
                         parents=fws[1],
                         gaussian_input_params=esp_gaussian_inputs,
                         working_dir=working_dir,
                         cart_coords=cart_coords,
                         spec={"proceed": {"has_gaussian_completed": True,
                                           "stationary_type": "Minimum"}},
                         **kwargs
                         )
    )
    return Workflow(fws,
                    name=get_job_name(mol, name),
                    **kwargs)


def get_nmr_tensors(mol_file=None,
                    smiles=None,
                    db=None,
                    name="nmr_tensor_calculation",
                    working_dir=None,
                    opt_gaussian_inputs=None,
                    freq_gaussian_inputs=None,
                    nmr_gaussian_inputs=None,
                    cart_coords=True,
                    save_to_db=False,
                    update_duplicates=False,
                    save_mol_file=False,
                    oxidation_states=None,
                    **kwargs):
    fws = []
    working_dir = working_dir or os.getcwd()
    mol, list_fws = common_fw(mol_file=mol_file,
                              smiles=smiles,
                              db=db,
                              working_dir=working_dir,
                              opt_gaussian_inputs=opt_gaussian_inputs,
                              freq_gaussian_inputs=freq_gaussian_inputs,
                              cart_coords=cart_coords,
                              save_to_db=save_to_db,
                              update_duplicates=update_duplicates,
                              oxidation_states=oxidation_states,
                              save_mol_file=save_mol_file,
                              **kwargs)
    fws += list_fws
    nmr_gaussian_inputs = nmr_gaussian_inputs or {}
    if "route_parameters" not in nmr_gaussian_inputs:
        nmr_gaussian_inputs.update({"route_parameters": {"NMR": "GIAO"}})
    fws.append(
        CalcFromRunsDBFW(db,
                         input_file="mol_nmr.com",
                         output_file="mol_nmr.out",
                         name=get_job_name(mol, "nmr"),
                         parents=fws[1],
                         gaussian_input_params=nmr_gaussian_inputs,
                         working_dir=working_dir,
                         cart_coords=cart_coords,
                         spec={"proceed": {"has_gaussian_completed": True,
                                           "stationary_type": "Minimum"}},
                         **kwargs
                         )
    )
    return Workflow(fws,
                    name=get_job_name(mol, name),
                    **kwargs)


def get_binding_energies():
    pass
