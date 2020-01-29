import os

from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianOutput
from infrastructure.gaussian.database import GaussianCalcDb


@explicit_serialize
class ConvertToMoleculeObject(FiretaskBase):
    """
    Reads a molecule from a file, converts it to a mol object,
    and saves it as dict to mongodb.
    Supported file formats include
        xyz|pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv,
        gaussian input (gjf|g03|g09|com|inp),
        Gaussian output (.out), and
        pymatgen's JSON serialized molecules.
    Requires openbabel to be installed.
    """
    optional_params = ['working_dir', 'mol_file', 'db']

    def run_task(self, fw_spec):
        # TODO: Give the option to not save the molecule to the database
        # TODO: Raise a warning that the molecule exists in the database
        working_dir = fw_spec. \
            get('working_dir', self.get('working_dir', os.getcwd()))
        file_name = fw_spec.get("mol_file", self.get('mol_file'))
        file_path = os.path.join(working_dir, file_name)
        mol = Molecule.from_file(file_path)
        mol_db = GaussianCalcDb(**fw_spec['db'])
        mol_db.insert_molecule(mol, update_duplicates=False)
        fw_spec['prev_calc_molecule'] = mol
        mol_db = GaussianCalcDb(**fw_spec['db'])
        fw_spec['run'] = {'smiles': mol_db.get_smiles(mol)}
        mol_db.insert_molecule(mol, update_duplicates=False)
