import os

from fireworks import Firework
from infrastructure.gaussian.firetasks.geo_transformation import \
    ConvertToMoleculeObject
from infrastructure.gaussian.firetasks.write_inputs import WriteInput
from infrastructure.gaussian.firetasks.run_calc import RunGaussianDirect
from infrastructure.gaussian.firetasks.parse_outputs import ProcessRun


class MolFileToInput(Firework):
    def __init__(self, mol_file, db=None, name='molToInput', parents=None,
                 working_dir=None, input_file='mol.com', save_to_db=False,
                 gaussian_input_params=None, cart_coords=None,
                 oxidation_states=None, **kwargs):
        t = []

        working_dir = working_dir or os.getcwd()
        t.append(
            ConvertToMoleculeObject(mol_file=mol_file, working_dir=working_dir,
                                    db=db, save_to_db=save_to_db)
        )

        t.append(
            WriteInput(input_file=input_file, working_dir=working_dir,
                       gaussian_input_params=gaussian_input_params,
                       cart_coords=cart_coords,
                       oxidation_states=oxidation_states)
        )
        super(MolFileToInput, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class RunCalc(Firework):
    def __init__(self, name='runCalc', parents=None, working_dir=None,
                 input_file='mol.com', output_file='mol.out', **kwargs):
        t = [RunGaussianDirect(working_dir=working_dir, input_file=input_file,
                               output_file=output_file)]

        super(RunCalc, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class SaveToDB(Firework):
    def __init__(self, db, name='saveToDb', parents=None, working_dir=None,
                 input_file='mol.com', output_file='mol.out',
                 fw_spec_field=None, **kwargs):
        t = [ProcessRun(db=db, working_dir=working_dir, input_file=input_file,
                        output_file=output_file, fw_spec_field=fw_spec_field)]

        super(SaveToDB, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)
