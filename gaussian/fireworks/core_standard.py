import os
import logging
from fireworks import Firework
from infrastructure.gaussian.firetasks.geo_transformation import \
    ConvertToMoleculeObject, RetrieveMoleculeObject, RetrieveGaussianOutput
from infrastructure.gaussian.firetasks.write_inputs import WriteInput
from infrastructure.gaussian.firetasks.run_calc import RunGaussianDirect
from infrastructure.gaussian.firetasks.parse_outputs import GaussianToDB

logger = logging.getLogger(__name__)


def common_tasks(db,
                 input_file,
                 output_file,
                 working_dir,
                 fw_spec_field,
                 gaussian_input_params,
                 cart_coords,
                 oxidation_states):
    return \
        [WriteInput(input_file=input_file,
                    working_dir=working_dir,
                    gaussian_input_params=gaussian_input_params,
                    cart_coords=cart_coords,
                    oxidation_states=oxidation_states),
         RunGaussianDirect(working_dir=working_dir,
                           input_file=input_file,
                           output_file=output_file),
         GaussianToDB(db=db,
                      working_dir=working_dir,
                      input_file=input_file,
                      output_file=output_file,
                      fw_spec_field=fw_spec_field)]


class CalcFromMolFileFW(Firework):
    def __init__(self,
                 mol_file,
                 db,
                 name="calc_from_mol_file",
                 parents=None,
                 working_dir=None,
                 input_file="mol.com",
                 output_file="mol.out",
                 gaussian_input_params=None,
                 cart_coords=True,
                 save_to_db=False,
                 update_duplicates=False,
                 oxidation_states=None,
                 fw_spec_field=None,
                 **kwargs):
        t = []
        working_dir = working_dir or os.getcwd()

        t.append(
            ConvertToMoleculeObject(db=db,
                                    mol_file=mol_file,
                                    working_dir=working_dir,
                                    save_to_db=save_to_db,
                                    update_duplicates=update_duplicates)
        )
        t += common_tasks(db,
                          input_file,
                          output_file,
                          working_dir,
                          fw_spec_field,
                          gaussian_input_params,
                          cart_coords,
                          oxidation_states)
        super(CalcFromMolFileFW, self).__init__(t,
                                                parents=parents,
                                                name=name,
                                                **kwargs)


class CalcFromRunsDBFW(Firework):
    def __init__(self,
                 db,
                 name="calc_from_runs_db",
                 parents=None,
                 gaussian_input_params=None,
                 working_dir=None,
                 input_file="mol.com",
                 output_file="mol.out",
                 cart_coords=True,
                 fw_spec_field=None,
                 **kwargs):
        t = []
        working_dir = working_dir or os.getcwd()

        t.append(
            RetrieveGaussianOutput(db=db,
                                   gaussian_input_params=gaussian_input_params)
        )
        t += common_tasks(db,
                          input_file,
                          output_file,
                          working_dir,
                          fw_spec_field,
                          gaussian_input_params,
                          cart_coords,
                          oxidation_states=None)
        super(CalcFromRunsDBFW, self).__init__(t,
                                               parents=parents,
                                               name=name,
                                               **kwargs)


class CalcFromMolDBFW(Firework):
    def __init__(self,
                 db,
                 smiles=None,
                 name="calc_from_mol_db",
                 parents=None,
                 gaussian_input_params=None,
                 working_dir=None,
                 input_file="mol.com",
                 output_file="mol.out",
                 cart_coords=True,
                 oxidation_states=None,
                 fw_spec_field=None,
                 **kwargs):
        t = []
        working_dir = working_dir or os.getcwd()

        t.append(
            RetrieveMoleculeObject(db=db,
                                   smiles=smiles)
        )
        t += common_tasks(db,
                          input_file,
                          output_file,
                          working_dir,
                          fw_spec_field,
                          gaussian_input_params,
                          cart_coords,
                          oxidation_states)
        super(CalcFromMolDBFW, self).__init__(t,
                                              parents=parents,
                                              name=name,
                                              **kwargs)
