from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.core.structure import Molecule
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
    # _fw_name = "Convert To Molecules Task"

    def run_task(self, fw_spec):
        # TODO: FIGURE OUT HOW TO SAVE TO THE DATABASE
        file_name = fw_spec["mol_file"]
        mol = Molecule.from_file(file_name)
        fw_spec['prev_calc_molecule'] = mol
        mol_db = GaussianCalcDb(**fw_spec['db'])
        fw_spec['run'] = {'smiles': mol_db.get_smiles(mol)}
        mol_db.insert_molecule(mol, update_duplicates=False)
