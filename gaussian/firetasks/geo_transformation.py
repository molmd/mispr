import os
import logging

from bson.objectid import ObjectId

from fireworks.core.firework import FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianInput

from infrastructure.gaussian.database import GaussianCalcDb
from infrastructure.gaussian.utils.utils import get_db, get_mol_from_file, \
    get_mol_from_db

logger = logging.getLogger(__name__)


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
        mol = get_mol_from_file(file_name, working_dir)
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
        mol_db = get_db(self.get('db'))
        mol = get_mol_from_db(smiles, mol_db)
        if mol and self.get("save_mol_file", False):
            working_dir = self.get('working_dir', os.getcwd())
            file_name = self.get(
                'filename.{}'.format(self.get('fmt', 'xyz')),
                'mol.{}'.format(self.get('fmt', 'xyz')))
            mol_file = os.path.join(working_dir, file_name)
            mol.to(self.get('fmt', 'xyz'), mol_file)
        fw_spec['prev_calc_molecule'] = mol


@explicit_serialize
class RetrieveGaussianOutput(FiretaskBase):
    """
    Returns a Gaussian output object from the database and converts it to a
    Gaussian input object
    """
    required_params = []
    optional_params = ["db", "gaussian_input_params", "smiles", "functional",
                       "basis", "type", "phase", "tag"]

    def run_task(self, fw_spec):
        # TODO: use alphabetical formula as a search criteria

        # get db credentials
        run_db = get_db(self.get('db'))

        # if a Gaussian output dict is passed through fw_spec
        if fw_spec.get("gaussian_output_id"):
            run = run_db.runs. \
                find_one({"_id": ObjectId(fw_spec.get("gaussian_output_id"))})
            proceed_keys = fw_spec.get("proceed", {})
            for k, v in proceed_keys.items():
                if run["output"].get(k, run["output"]["output"].get(k)) != v:
                    raise ValueError(
                        f"The condition for {k} is not met, Terminating"
                    )

        # if a Gaussian output dictionary is retrieved from db
        else:
            smiles = self['smiles']
            functional = self.get('functional')
            basis = self.get('basis')
            type_ = self.get('type')
            phase = self.get('phase')
            if 'tag' in self:
                tag = self['tag']
                run = run_db.retrieve_run(smiles, type_, functional, basis,
                                          phase, tag=tag)
            else:
                run = run_db.retrieve_run(smiles, type_, functional, basis,
                                          phase)
            if not run:
                raise Exception("Gaussian output is not in the database")
            run = max(run, key=lambda i: i['last_updated'])

        # create a gaussian input object from run
        if self.get("gaussian_input_params") is None:
            logger.info("No gaussian input parameters provided; will use "
                        "run parameters")
        inputs = {}
        for k, v in run['input'].items():
            # use gaussian_input_params if defined, otherwise use run parameters
            inputs[f'{k}'] = self.get("gaussian_input_params", {}).\
                get(f'{k}', run['input'].get(f'{k}'))
        inputs['molecule'] = run['output']['output']['molecule']
        gaussin = GaussianInput.from_dict(inputs)
        fw_spec["gaussian_input"] = gaussin


@explicit_serialize
class AttachFunctionalGroup(FiretaskBase):
    """
    Attaches a functional group to a molecule; requires the name of the functional
    group and the smiles representation of the molecule to be read from the database
    (can be taken from both the molecule collection or the runs collection)
    """
    required_params = ["func_grp", "index"]
    optional_params = ["db", "smiles", "bond_order", "save_to_db",
                       "update_duplicates",
                       "save_mol_file", "fmt", "filename", "working_dir"]
    # TODO: check if it is better to split this into multiple firetasks (one of
    # which has already been created above (Retrieve Molecule from db)

    def run_task(self, fw_spec):
        db = get_db(self.get('db'))
        if fw_spec.get("prev_calc_molecule"):
            mol = fw_spec.get("prev_calc_molecule")
        elif self.get("smiles"):
            mol_dict = db.retrieve_molecule(self.get["smiles"])
            mol = Molecule.from_dict(mol_dict)
        func_grp = self['func_grp']
        func_grp_dict = db.retrieve_fg(func_grp)
        func_group = Molecule(func_grp_dict["species"], func_grp_dict["coords"])
        derived_mol = Molecule.copy(mol)
        derived_mol.substitute(index=self['index'], func_grp=func_group)
        if self.get('save_to_db', True):
            db.insert_derived_mol(derived_mol,
                                  update_duplicates=self.get(
                                      'update_duplicates', False))
        if self.get("save_mol_file", False):
            working_dir = self.get('working_dir', os.getcwd())
            file_name = self.get(
                'filename.{}'.format(self.get('fmt', 'xyz')),
                'mol.{}'.format(self.get('fmt', 'xyz')))
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
        # TODO: take common things from functions
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





