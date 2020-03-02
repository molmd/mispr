import os
from pymatgen.core.structure import Molecule
from pymatgen.io.babel import BabelMolAdaptor
from fireworks.utilities.fw_utilities import get_slug
from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks import FileWriteTask
import pybel as pb


def get_mol_from_file(mol_file, working_dir=None):
    working_dir = working_dir or os.getcwd()
    file_path = os.path.join(working_dir, mol_file)
    mol = Molecule.from_file(file_path)
    return mol


def get_mol_from_db(smiles, db):
    mol_db = get_db(db)
    mol_dict = mol_db.retrieve_molecule(smiles)
    if mol_dict is None:
        raise Exception("Molecule is not found in the database")
    mol = Molecule.from_dict(mol_dict)
    return mol


def get_chem_schema(mol):
    mol_dict = mol.as_dict()
    comp = mol.composition
    a = BabelMolAdaptor(mol)
    pm = pb.Molecule(a.openbabel_mol)
    # svg = pm.write("svg")
    mol_dict.update({"smiles": pm.write("can").strip(),
                     "formula": comp.formula,
                     "formula_pretty": comp.reduced_formula,
                     "formula_anonymous": comp.anonymized_formula,
                     "formula_alphabetical": comp.alphabetical_formula,
                     "chemsys": comp.chemical_system,
                     "nsites": mol.num_sites,
                     "nelements":
                         len(comp.chemical_system.replace('-', ' ').split(' ')),
                     "is_ordered": mol.is_ordered,
                     "is_valid": mol.is_valid()})
    return mol_dict


def get_mol_formula(mol):
    mol_schema = get_chem_schema(mol)
    return mol_schema["formula_pretty"]


def get_job_name(mol, name):
    mol_formula = get_mol_formula(mol)
    job_name = "{}_{}".format(mol_formula, name)
    return job_name


def add_namefile(original_wf, use_slug=True):
    for idx, fw in enumerate(original_wf.fws):
        fname = "FW--{}".format(fw.name)
        if use_slug:
            fname = get_slug(fname)

        t = FileWriteTask(files_to_write=[{"filename": fname, "contents": ""}])
        original_wf.fws[idx].tasks.insert(0, t)
    return original_wf


def get_db(input_db=None):
    from infrastructure.gaussian.database import GaussianCalcDb
    if not input_db:
        input_db = f"{CONFIG_FILE_DIR}/db.json"
        if not os.path.isfile(input_db):
            raise FileNotFoundError("Please provide the database configurations")
    if isinstance(input_db, dict):
        db = GaussianCalcDb(**input_db)
    else:
        db = GaussianCalcDb.from_db_file(input_db)

    return db
