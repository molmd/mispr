import os
import logging

from bson.objectid import ObjectId

from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.core.structure import Molecule

from infrastructure.gaussian.utils.utils import get_db, process_mol, \
    get_run_from_fw_spec, process_run

logger = logging.getLogger(__name__)

DEFAULT_KEY = 'gout_key'


@explicit_serialize
class ProcessMoleculeInput(FiretaskBase):
    required_params = ["mol"]
    optional_params = ["operation_type", "db", "working_dir", "save_to_db",
                       "update_duplicates", "save_to_file", "fmt", "filename",
                       "from_fw_spec"]

    @staticmethod
    def _run_to_mol_object(run):
        return Molecule.from_dict(run['output']['output']['molecule'])

    @staticmethod
    def _from_fw_spec(mol, fw_spec):
        available_runs = fw_spec['gaussian_output']
        if not isinstance(mol, dict):
            mol = ProcessMoleculeInput._run_to_mol_object(available_runs[mol])
        else:
            if isinstance(mol['mol'], list):
                mol['mol'] = \
                    [ProcessMoleculeInput._run_to_mol_object(available_runs[i])
                     for i in mol['mol']]
            else:
                mol['mol'] = ProcessMoleculeInput.\
                    _run_to_mol_object(available_runs[mol['mol']])
        return mol

    def run_task(self, fw_spec):
        mol = self["mol"]
        operation_type = self.get("operation_type", "get_from_mol")
        working_dir = self.get('working_dir', os.getcwd())
        db = self.get('db')

        if self.get('from_fw_spec'):
            mol = self._from_fw_spec(mol, fw_spec)

        output_mol = process_mol(operation_type=operation_type, mol=mol,
                                 working_dir=working_dir, db=db)

        if self.get("save_to_db"):
            db = get_db(db) if db else get_db()
            update_duplicates = self.get("update_duplicates", False)
            db.insert_molecule(output_mol, update_duplicates=update_duplicates)

        if self.get("save_to_file"):
            fmt = self.get("fmt", 'xyz')
            filename = self.get("filename", 'mol')
            file = os.path.join(working_dir, f"{filename}.{fmt}")
            output_mol.to(fmt, file)
        fw_spec['prev_calc_molecule'] = output_mol # Note: This should ideally
        # be part of FWaction, however because mpi doesn't support pymatgen, we
        # should be careful about what is being passed to the next firework


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
        mol = process_mol(file_name, working_dir)
        if self.get('save_to_db', True):
            mol_db = get_db(self.get('db'))
            mol_db.insert_molecule(
                mol,
                update_duplicates=self.get('update_duplicates', False)
            )
        fw_spec['prev_calc_molecule'] = mol  # Note: This should ideally be
        # part of FWaction, however because mpi doesn't support pymatgen, we
        # should be careful about what is being passed to the next firework


@explicit_serialize
class RetrieveMoleculeObject(FiretaskBase):
    """
    Returns a molecule object from the database using the smiles as an identifier
    """
    required_params = ["smiles"]
    optional_params = ["db", "save_mol_file", "fmt", "filename", "working_dir"]

    def run_task(self, fw_spec):
        # TODO: use alphabetical formula as a search criteria
        smiles = self['smiles']
        mol = process_mol(smiles, self.get('db'))
        if mol and self.get("save_mol_file", False):
            working_dir = self.get('working_dir', os.getcwd())
            file_name = self.get(
                'filename.{}'.format(self.get('fmt', 'xyz')),
                'mol.{}'.format(self.get('fmt', 'xyz')))
            mol_file = os.path.join(working_dir, file_name)
            mol.to(self.get('fmt', 'xyz'), mol_file)
        fw_spec['prev_calc_molecule'] = mol


@explicit_serialize
class AttachFunctionalGroup(FiretaskBase):
    """
    Attaches a functional group to a molecule; requires the name of the functional
    group and the smiles representation of the molecule to be read from the database
    (can be taken from both the molecule collection or the runs collection)
    """
    required_params = ["func_grp", "index"]
    optional_params = ["db", "molecule", "bond_order", "save_to_db",
                       "update_duplicates",
                       "save_mol_file", "fmt", "filename", "working_dir"]
    # TODO: check if it is better to split this into multiple firetasks (one of
    # which has already been created above (Retrieve Molecule from db)

    def run_task(self, fw_spec):
        db = get_db(self.get('db'))
        if fw_spec.get("prev_calc_molecule"):
            mol = fw_spec.get("prev_calc_molecule")
        elif self.get("molecule"):
            mol = self.get("molecule")
        func_grp_dict = db.retrieve_fg(self['func_grp'])
        func_grp = Molecule(func_grp_dict["species"], func_grp_dict["coords"])
        derived_mol = Molecule.copy(mol)
        derived_mol.substitute(index=self['index'], func_grp=func_grp)
        if self.get('save_to_db', True):
            db.insert_derived_mol(derived_mol,
                                  update_duplicates=self.get(
                                      'update_duplicates', False))
        if self.get("save_mol_file", False):
            working_dir = self.get('working_dir', os.getcwd())
            file_name = self.get('filename', 'derived_mol')
            file_name = '{}.{}'.format(file_name, self.get('fmt', 'xyz')),
            derived_mol_file = os.path.join(working_dir, file_name)
            derived_mol.to(self.get('fmt', 'xyz'), derived_mol_file)
        fw_spec['prev_calc_molecule'] = derived_mol


@explicit_serialize
class LinkMolecules(FiretaskBase):
    """
    Links two molecules using one site from the first and another site from the
    second molecule. Currently takes the molecules from the db using their
    smiles representation.
    """
    required_params =["index1", "index2"]
    optional_params = ["db", "smiles1", "smiles2", "bond_order", "save_to_db",
                       "update_duplicates", "save_mol_file", "fmt",
                       "filename", "working_dir"]

    def run_task(self, fw_spec):
        # TODO: take mol1 and mol2 from previous calculations
        db = get_db(self.get("db"))
        mol1_dict = db.retrieve_molecule(self.get["smiles1"])
        mol1 = Molecule.from_dict(mol1_dict)
        mol2_dict = db.retrieve_molecule(self.get["smiles2"])
        mol2 = Molecule.from_dict(mol2_dict)
        linked_mol = mol1.link(mol2, self["index1"], self["index2"],
                               self.get["bond_order"])
        if self.get("save_to_db", True):
            db.insert_molecule(linked_mol, update_duplicates=self.
                               get("update_duplicates", False))
        if self.get("save_mol_file", False):
            working_dir = self.get("working_dir", os.getcwd())
            file_name = self.get("filename.{}".format(self.get("fmt", "xyz")),
                                 "mol.{}".format(self.get("fmt", "xyz")))
            linked_mol_file = os.path.join(working_dir, file_name)
            linked_mol.to(self.get("fmt", "xyz"), linked_mol_file)
        fw_spec["prev_calc_molecule"] = linked_mol
