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
    required_params = ['mol_file']
    optional_params = ['db', 'working_dir', 'save_to_db', 'update_duplicates']

    def run_task(self, fw_spec):
        working_dir = self.get('working_dir', os.getcwd())
        file_name = self["mol_file"]
        file_path = os.path.join(working_dir, file_name)
        mol = Molecule.from_file(file_path)
        if self.get('save_to_db', True):
            mol_db = GaussianCalcDb(**self.get('db'))
            mol_db.insert_molecule(mol, update_duplicates=self.
                                   get('update_duplicates', False))
        fw_spec['prev_calc_molecule'] = mol  # Note: This should ideally be
        # part of FWaction, however because mpi doesn't support pymatgen, we
        # should be careful about what is being passed to the next firework


@explicit_serialize
class RetrieveMoleculeObject(FiretaskBase):
    """
    Returns a molecule object from the database using the smiles as an identifier
    """
    optional_params = ['smiles', 'db']

    def run_task(self, fw_spec):
        smiles = fw_spec.get('smiles')
        mol_db = GaussianCalcDb(**fw_spec['db'])
        mol = mol_db.retrieve_molecule(smiles)
        if mol is None:
            raise Exception("Molecule is not found in the database")
        fw_spec['prev_calc_molecule'] = mol


@explicit_serialize
class RetrieveGaussianOutput(FiretaskBase):
    """
    Returns a Gaussian output object from the database and converts it to a
    Gaussian input object
    """
    optional_params = ['smiles', 'functional', 'basis', 'type', 'db',
                       "gaussian_input_params"]

    def run_task(self, fw_spec):
        smiles = fw_spec.get('smiles')
        functional = fw_spec.get('functional', 'B3LYP')
        basis = fw_spec.get('basis', '6-31+G*')
        type_ = fw_spec.get('type', 'Opt')
        run_db = GaussianCalcDb(**fw_spec['db'])
        if 'tag' in fw_spec:
            tag = fw_spec['tag']
            run = run_db.retrieve_run(smiles, type_, functional, basis, tag=tag)
        else:
            run = run_db.retrieve_run(smiles, type_, functional, basis)
        if not run:
            raise Exception("Gaussian output is not in the database")
        run = max(run, key=lambda i: i['last_updated'])

        run['output']['charge'] = self.get("gaussian_input_params", {}).\
            get('charge', run['output'].get('charge'))
        run['output']['spin_multiplicity'] = \
            self.get("gaussian_input_params", {}).\
                get('spin_multiplicity', run['output'].get('spin_multiplicity'))
        run['output']['title'] = self.get("gaussian_input_params", {}).\
            get('title', run['output'].get('title'))
        run['output']['input'] = {**run['output'].get('input', {}),
                                  **self.get('gaussian_input_params', {})}
        gaussin = GaussianOutput.from_dict_to_input(run['output'])
        fw_spec["gaussian_input"] = gaussin


@explicit_serialize
class AttachFunctionalGroup(FiretaskBase):
    """
    Attaches a functional group to a molecule; requires the name of the functional
    group and the smiles representation of the molecule to be read from the database
    (can be taken from both the molecule collection or the runs collection)
    """

    def run_task(self, fw_spec):
        pass
