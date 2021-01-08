import os
import logging

from fireworks import Firework, Workflow

from infrastructure.gaussian.utils.utils import process_mol, get_job_name, \
    get_mol_formula, process_run
from infrastructure.gaussian.firetasks.parse_outputs import ProcessRun
from infrastructure.gaussian.fireworks.core_standard import CalcFromMolFW, \
    CalcFromRunsDBFW

logger = logging.getLogger(__name__)

WORKFLOW_KWARGS = Workflow.__init__.__code__.co_varnames


# TODO: avoid overwriting a directory if it exists (happens when molecules have
#  same mol formula)
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
              check_result=None,
              **kwargs):
    # TODO: checking criteria not working
    #   TODO: mol is not returned if process_mol_func is set to False
    #    TODO: this function is changed to return three outputs, should
    #    change everything

    def _recursive_key_check(grun, grun_key):
        key_exists = grun_key in grun
        if not key_exists:
            for key, value in grun.items():
                if isinstance(value, dict):
                    key_exists = key_exists or _recursive_key_check(value,
                                                                    grun_key)
        return key_exists

    fws = []

    if not gout_key:
        gout_key = "mol"
    if not skip_opt_freq:
        if process_mol_func:
            mol = process_mol(mol_operation_type, mol, db=db,
                              working_dir=working_dir)
            mol_operation_type = "get_from_mol"
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
        dir_struct = [dir_head] + kwargs.get("dir_structure", [])
        working_dir = os.path.join(working_dir, *dir_struct)
        opt_input_file = f"{label}_opt.com"
        opt_output_file = f"{label}_opt.out"
        freq_input_file = f"{label}_freq.com"
        freq_output_file = f"{label}_freq.out"

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
                               gout_key=gout_key + "_opt",
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
            # TODO: add option for user to add label or mol_name like above
            run = process_run(mol_operation_type, mol, db=db,
                              working_dir=working_dir)
            if check_result:
                keys_exist = []
                for grun_key in check_result:
                    keys_exist.append(_recursive_key_check(run, grun_key))
                not_found_ind = [ind for ind, boolean in enumerate(keys_exist)
                                 if not boolean]
                if not_found_ind:
                    not_found_keys = [check_result[ind] for ind in not_found_ind]
                    raise ValueError("Gaussian output does not include {}. "
                                     "Stopping.".format(not_found_keys))
            mol = process_mol("get_from_run_dict", run, db=db,
                              working_dir=working_dir)
            spec = kwargs.pop("spec", {})
            label = get_mol_formula(mol)
            if not dir_head:
                dir_head = label
            dir_struct = [dir_head] + kwargs.get("dir_structure", [])
            working_dir = os.path.join(working_dir, *dir_struct)
            if "tag" in kwargs:
                spec.update({"tag": kwargs["tag"]})
            spec.update({"_launch_dir": working_dir})
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
    return mol, label, fws
