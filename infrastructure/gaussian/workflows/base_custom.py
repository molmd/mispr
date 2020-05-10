from fireworks import Workflow

from infrastructure.gaussian.fireworks.core_custom import MolFileToInput, \
    RunCalc, SaveToDB


def get_esp_charges(mol_file, db, working_dir=None,
                                input_file='mol.com', output_file='mol.out',
                                gaussian_input_params=None, cart_coords=None,
                                save_to_db=False, oxidation_states=None,
                                **kwargs):
    firework1 = MolFileToInput(mol_file, db, working_dir=working_dir,
                               input_file=input_file, save_to_db=save_to_db,
                               gaussian_input_params=gaussian_input_params,
                               cart_coords=cart_coords,
                               oxidation_states=oxidation_states, **kwargs)

    spec = kwargs.pop('spec', {})
    spec.update({'_queueadapter': {'queue': 'mpi'}})
    firework2 = RunCalc(parents=[firework1], working_dir=working_dir,
                        input_file=input_file, output_file=output_file,
                        spec=spec, **kwargs)

    spec.update({'_queueadapter': {'queue': 'rajputlab'}})
    firework3 = SaveToDB(db, parents=[firework2], working_dir=working_dir,
                         input_file=input_file, output_file=output_file,
                         spec=spec, **kwargs)

    return Workflow([firework1, firework2, firework3])