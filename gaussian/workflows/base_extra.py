import os
import logging

from fireworks import Firework, Workflow
from pymatgen.core.structure import Molecule
from infrastructure.gaussian.firetasks.parse_outputs import \
    BindingEnergytoDB
from bson.objectid import ObjectId
from infrastructure.gaussian.firetasks.parse_outputs import ProcessRun
from infrastructure.gaussian.fireworks.core_standard import CalcFromMolFW, \
    CalcFromRunsDBFW
from infrastructure.gaussian.utils.utils import process_mol, get_job_name, \
    get_mol_formula, process_run, get_db, get_run_from_fw_spec
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize

logger = logging.getLogger(__name__)

WORKFLOW_KWARGS = Workflow.__init__.__code__.co_varnames

STANDARD_OPT_GUASSIAN_INPUT = {"functional": "B3LYP",
                               "basis_set": "6-31G(d)",
                               "route_parameters": {"Opt": None},
                               "link0_parameters": {"%chk": "checkpoint.chk",
                                                    "%mem": "1000MW",
                                                    "%NProcShared": "12"}}

MOL_LIST_ORDER = ["pristine", "reduced", "oxidized"]

@explicit_serialize
class RedoxtoDB(FiretaskBase):
    required_params = ["keys", 'main_run_key', 'new_prop']
    optional_params = ["db"]

    def run_task(self, fw_spec):
        run_db = get_db(self.get('db'))
        runs = [get_run_from_fw_spec(fw_spec, i, run_db) for i in self['keys']]
        props = [i['output']['output']["final_energy"] for i in runs]
        result = (props[2] - (props[0] + props[1])) * 27.2114
        main_run = get_run_from_fw_spec(fw_spec, self['main_run_key'], run_db)
        run_db.update_run(new_values={self["new_prop"]: result},
                          _id=ObjectId(main_run['_id']))
        logger.info("binding energy calculation complete")


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
            if "tag" in kwargs:
                spec.update({'tag': kwargs["tag"]})
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
    return mol, fws


def ss_add_h_redox(mol_operation_type,
                   mol_list,
                   db=None,
                   name="ss_add_h_redox_calculation",
                   working_dir=None,
                   opt_gaussian_inputs=None,
                   freq_gaussian_inputs=None,
                   cart_coords=True,
                   oxidation_states=None,
                   skip_opt_freq=False,
                   **kwargs):
    fws = []
    working_dir = working_dir or os.getcwd()
    for index, mol in mol_list:
        cur_working_dir = os.path.join(working_dir, MOL_LIST_ORDER[index])
        os.makedirs(cur_working_dir, exist_ok=True)
        mol, opt_freq_fws = common_fw(mol_operation_type=mol_operation_type,
                                      mol=mol,
                                      working_dir=cur_working_dir,
                                      db=db,
                                      opt_gaussian_inputs=opt_gaussian_inputs,
                                      freq_gaussian_inputs=freq_gaussian_inputs,
                                      cart_coords=cart_coords,
                                      oxidation_states=oxidation_states,
                                      skip_opt_freq=skip_opt_freq,
                                      **kwargs)
        fws += opt_freq_fws

    return Workflow(fws,
                    name=get_job_name(mol, name),
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})