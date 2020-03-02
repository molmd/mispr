import os
import logging

from fireworks import Workflow
from infrastructure.gaussian.fireworks.core_standard import CalcFromMolFileFW, \
    CalcFromRunsDBFW, CalcFromMolDBFW
from infrastructure.gaussian.utils.utils import get_mol_from_file, \
    get_mol_from_db, get_job_name

logger = logging.getLogger(__name__)

STANDARD_OPT_GUASSIAN_INPUT = {"functional": "B3LYP",
                               "basis_set": "6-31G(d)",
                               "route_parameters": {"Opt": None},
                               "link0_parameters": {"%chk": "checkpoint.chk",
                                                    "%mem": "1000MW",
                                                    "%NProcShared": "24"}}


def common_fw(working_dir,
              opt_gaussian_inputs,
              freq_gaussian_inputs,
              cart_coords,
              save_to_db,
              update_duplicates,
              oxidation_states,
              db=None,
              mol_file=None,
              smiles=None,
              **kwargs):
    opt_gaussian_inputs = opt_gaussian_inputs or {}
    opt_gaussian_inputs = {**STANDARD_OPT_GUASSIAN_INPUT, **opt_gaussian_inputs}
    if "opt" not in [i.lower() for i in
                     opt_gaussian_inputs["route_parameters"]]:
        raise ValueError("The Opt keyword is missing from the input file")
    if mol_file is not None:
        mol = get_mol_from_file(mol_file, working_dir)
        fw1 = CalcFromMolFileFW(mol_file,
                                db,
                                name=get_job_name(mol, "optimization"),
                                working_dir=working_dir,
                                input_file="mol_opt.com",
                                output_file="mol_opt.out",
                                gaussian_input_params=opt_gaussian_inputs,
                                cart_coords=cart_coords,
                                save_to_db=save_to_db,
                                update_duplicates=update_duplicates,
                                oxidation_states=oxidation_states,
                                **kwargs
                                )
    elif smiles is not None:
        mol = get_mol_from_db(smiles, db)
        fw1 = CalcFromMolDBFW(db,
                              smiles,
                              name=get_job_name(mol, "optimization"),
                              working_dir=working_dir,
                              input_file="mol_opt.com",
                              output_file="mol_opt.out",
                              gaussian_input_params=opt_gaussian_inputs,
                              cart_coords=cart_coords,
                              oxidation_states=oxidation_states,
                              **kwargs
                              )
    else:
        raise ValueError("No molecule provided. Either a molecule file or a "
                         "smiles representation from the molecules collection "
                         "should be provided as an input")
    # if no freq_gaussian_inputs are provided, parameters from prev opt are used
    # except for the route parameters which are replaced with the Freq keyword
    freq_gaussian_inputs = freq_gaussian_inputs or {}
    if "route_parameters" not in freq_gaussian_inputs:
        freq_gaussian_inputs.update({"route_parameters": {"Freq": None}})
    if "freq" not in [i.lower() for i in
                      freq_gaussian_inputs["route_parameters"]]:
        raise ValueError('The Freq keyword is missing from the input file')
    fw2 = CalcFromRunsDBFW(db,
                           input_file="mol_freq.com",
                           output_file="mol_freq.out",
                           name=get_job_name(mol, "frequency"),
                           parents=fw1,
                           gaussian_input_params=freq_gaussian_inputs,
                           working_dir=working_dir,
                           cart_coords=cart_coords,
                           spec={"proceed": {"has_gaussian_completed": True}},
                           **kwargs
                           )
    return mol, [fw1, fw2]


def get_esp_charges(mol_file=None,
                    smiles=None,
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
                    oxidation_states=None,
                    **kwargs):
    fws = []
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
