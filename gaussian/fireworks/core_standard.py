import os
import logging
from fireworks import Firework
from infrastructure.gaussian.firetasks.geo_transformation import \
    ProcessMoleculeInput
from infrastructure.gaussian.firetasks.write_inputs import WriteInput
from infrastructure.gaussian.firetasks.run_calc import RunGaussianDirect
from infrastructure.gaussian.firetasks.parse_outputs import ProcessRun, \
    RetrieveGaussianOutput

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = ["spec", "name", "launches", "archived_launches", "state",
                   "created_on", "fw_id", "parents", "updated_on"]


def common_tasks(db,
                 input_file,
                 output_file,
                 working_dir,
                 gaussian_input_params,
                 cart_coords,
                 oxidation_states,
                 **kwargs):

    return \
        [WriteInput(input_file=input_file,
                    working_dir=working_dir,
                    gaussian_input_params=gaussian_input_params,
                    cart_coords=cart_coords,
                    oxidation_states=oxidation_states,
                    **{i: j for i, j in kwargs.items() if i in
                       WriteInput.required_params +
                       WriteInput.optional_params}),
         RunGaussianDirect(working_dir=working_dir,
                           input_file=input_file,
                           output_file=output_file,
                           **{i: j for i, j in kwargs.items() if i in
                              RunGaussianDirect.required_params +
                              RunGaussianDirect.optional_params}),
         ProcessRun(run=output_file,
                    operation_type="get_from_file",
                    db=db,
                    working_dir=working_dir,
                    input_file=input_file,
                    **{i: j for i, j in kwargs.items() if i in
                       ProcessRun.required_params +
                       ProcessRun.optional_params})]


class CalcFromMolFW(Firework):
    def __init__(self,
                 mol,
                 mol_operation_type="get_from_mol",
                 db=None,
                 name="calc_from_mol",
                 parents=None,
                 working_dir=None,
                 input_file="mol.com",
                 output_file="mol.out",
                 gaussian_input_params={},
                 cart_coords=True,
                 oxidation_states=None,
                 tag="unknown",
                 **kwargs):
        t = []
        working_dir = working_dir or os.getcwd()

        t.append(ProcessMoleculeInput(mol=mol,
                                      operation_type=mol_operation_type,
                                      db=db,
                                      working_dir=working_dir,
                                      **{i: j for i, j in kwargs.items() if i in
                                         ProcessMoleculeInput.required_params +
                                         ProcessMoleculeInput.optional_params}
                                      )
                 )

        t += common_tasks(db,
                          input_file,
                          output_file,
                          working_dir,
                          gaussian_input_params,
                          cart_coords,
                          oxidation_states,
                          **kwargs)
        spec = kwargs.pop('spec', {})
        spec.update({'tag': tag})
        super(CalcFromMolFW, self).__init__(t,
                                            parents=parents,
                                            name=name,
                                            spec=spec,
                                            **{i: j for i, j in kwargs.items()
                                               if i in FIREWORK_KWARGS})


class CalcFromRunsDBFW(Firework):
    def __init__(self,
                 db=None,
                 name="calc_from_runs_db",
                 parents=None,
                 gaussian_input_params=None,
                 working_dir=None,
                 input_file="mol.com",
                 output_file="mol.out",
                 cart_coords=True,
                 tag="unknown",
                 **kwargs):
        t = []
        working_dir = working_dir or os.getcwd()

        t.append(
            RetrieveGaussianOutput(db=db,
                                   gaussian_input_params=gaussian_input_params,
                                   **{i: j for i, j in kwargs.items() if i in
                                      RetrieveGaussianOutput.required_params +
                                      RetrieveGaussianOutput.optional_params}
                                   )
        )
        t += common_tasks(db,
                          input_file,
                          output_file,
                          working_dir,
                          gaussian_input_params,
                          cart_coords,
                          oxidation_states=None,
                          **kwargs)
        spec = kwargs.pop('spec', {})
        spec.update({'tag': tag})
        super(CalcFromRunsDBFW, self).__init__(t,
                                               parents=parents,
                                               name=name,
                                               spec=spec,
                                               **{i: j for i, j in
                                                  kwargs.items() if i in
                                                  FIREWORK_KWARGS})
