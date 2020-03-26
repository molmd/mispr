import os
import logging
import json
from bson.objectid import ObjectId
from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.babel import BabelMolAdaptor
from fireworks.utilities.fw_utilities import get_slug
from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks import FileWriteTask
import pybel as pb

logger = logging.getLogger(__name__)


def process_mol(operation_type, mol, **kwargs):
    working_dir = kwargs["working_dir"] if "working_dir" in kwargs \
        else os.getcwd()
    db = get_db(kwargs["db"]) if "db" in kwargs else get_db()

    if operation_type == 'get_from_mol':
        if not isinstance(mol, Molecule):
            raise Exception("mol is not a Molecule object; either "
                            "provide a Molecule object or use another "
                            "operation type with its corresponding inputs")
        output_mol = mol

    elif operation_type == 'get_from_file':
        # mol = file_path
        file_path = os.path.join(working_dir, mol)
        if not os.path.exists(file_path):
            raise Exception("mol is not a valid path; either provide a valid "
                            "path or use another operation type with its "
                            "corresponding inputs")

        output_mol = Molecule.from_file(file_path)

    elif operation_type == 'get_from_smiles':
        # mol = mol_smile
        mol_dict = db.retrieve_molecule(mol)
        if not mol_dict:
            raise Exception(
                "mol is not found in the database")
        output_mol = Molecule.from_dict(mol_dict)

    elif operation_type == 'get_from_gout':
        if not isinstance(mol, GaussianOutput):
            raise Exception("mol is not a GaussianOutput object; either "
                            "provide a GaussianOutput object or use another "
                            "operation type with its corresponding inputs")
        output_run = mol.as_dict()
        output_mol = Molecule.from_dict(output_run["output"]["molecule"])

    elif operation_type == 'get_from_run_dict':
        if not isinstance(mol, dict) and "output" not in mol:
            raise Exception("mol is not a GaussianOutput dictionary; either"
                            "provide a GaussianOutput dictionary or use another"
                            "operation type with its corresponding inputs")
        # TODO: sometimes it is mol['output']['output']['molecule']
        mol_dict = mol['output']['molecule']
        output_mol = Molecule.from_dict(mol_dict)

    elif operation_type == 'get_from_run_id':
        # mol = run_id
        run = db.runs.find_one({"_id": ObjectId(mol)})
        if not run:
            raise Exception("Gaussian run is not in the database")
        mol_dict = run['output']['output']['molecule']
        output_mol = Molecule.from_dict(mol_dict)

    elif operation_type == 'get_from_run_query':
        # mol = {'smiles': smiles, 'type': type, 'functional': func,
        #        'basis': basis, 'phase': phase, ...}
        logger.info("If the query criteria satisfy more than "
                    "one document, the last updated one will "
                    "be used. To perform a more specific "
                    "search, provide the document id using "
                    "gout_id")
        run = db.retrieve_run(**mol)
        if not run:
            raise Exception("Gaussian run is not in the database")
        run = max(run, key=lambda i: i['last_updated'])
        mol_dict = run['output']['output']['molecule']
        output_mol = Molecule.from_dict(mol_dict)

    elif operation_type == 'derive_molecule':
        # mol = {'operation_type': 'get_from_file', 'mol': file_path,
        #        'func_grp': func_group_name, ....}

        func_grp_name = mol.get('func_grp')
        if not func_grp_name:
            raise Exception("No FG provided; Provide the name of the FG "
                            "to be retrieved from the database")
        fg_dict = db.retrieve_fg(func_grp_name)
        if not fg_dict:
            raise Exception("FG is not found in the database")
        fg = Molecule(fg_dict["species"], fg_dict["coords"])

        output_mol = process_mol(operation_type=mol['operation_type'],
                                 mol=mol['mol'], **kwargs)
        output_mol.substitute(mol["index"], fg, mol["bond_order"])

    elif operation_type == 'link_molecules':
        # mol = {'operation_type': ['get_from_file', 'get_from_smiles'],
        #        'mol': ['kes/rasa/defe.xyz', 'mol_smiles'],
        #        'index': [3, 5],
        #        'bond_order': 2}

        # mol = {'operation_type': ['get_from_file', 'derive_molecule'],
        #        'mol': ['kes/rasa/defe.xyz', {'operation_type':
        #        'get_from_smiles, 'mol': smile}],
        #        'index': [3, 5],
        #        'bond_order': 2}
        linking_mol = process_mol(operation_type=mol['operation_type'][0],
                                  mol=mol['mol'][0], **kwargs)
        linked_mol = process_mol(operation_type=mol['operation_type'][1],
                                 mol=mol['mol'][1], **kwargs)
        output_mol = linking_mol.link(linked_mol, mol["index"][0],
                                      mol["index"][1], mol["bond_order"])
    else:
        raise ValueError(f'operation type {operation_type} is not supported')

    return output_mol


def process_run(operation_type, run, **kwargs):
    working_dir = kwargs["working_dir"] if "working_dir" in kwargs \
        else os.getcwd()
    db = get_db(kwargs["db"]) if "db" in kwargs else get_db()

    if operation_type == 'get_from_gout':
        if not isinstance(run, GaussianOutput):
            raise Exception("run is not a GaussianOutput object; either "
                            "provide a GaussianOutput object or use another "
                            "operation type with its corresponding inputs")
        output_run = run.as_dict()

    elif operation_type == 'get_from_file':
        # run = file_path
        file_path = os.path.join(working_dir, run)
        if not os.path.exists(file_path):
            raise Exception("run is not a valid path; either provide a valid "
                            "path or use another operation type with its "
                            "corresponding inputs")

        output_run = GaussianOutput(file_path).as_dict()

    elif operation_type == 'get_from_run_dict':
        if not isinstance(run, dict) and "output" not in run:
            raise Exception("run is not a GaussianOutput dictionary; either"
                            "provide a GaussianOutput dictionary or use another"
                            "operation type with its corresponding inputs")
        output_run = run

    elif operation_type == 'get_from_run_id':
        # run = run_id
        run = db.runs.find_one({"_id": ObjectId(run)})
        if not run:
            raise Exception("Gaussian run is not in the database")
        output_run = run['output']

    elif operation_type == 'get_from_run_query':
        # run = {'smiles': smiles, 'type': type, 'functional': func,
        #        'basis': basis, 'phase': phase, ...}
        logger.info("If the query criteria satisfy more than "
                    "one document, the last updated one will "
                    "be used. To perform a more specific "
                    "search, provide the document id using "
                    "gout_id")
        run = db.retrieve_run(**run)
        if not run:
            raise Exception("Gaussian run is not in the database")
        run = max(run, key=lambda i: i['last_updated'])
        output_run = run['output']

    else:
        raise ValueError(f'operation type {operation_type} is not supported')
    output_run = json.loads(json.dumps(output_run))
    output_run = recursive_signature_remove(output_run)
    return output_run


def get_run_from_fw_spec(fw_spec, key, run_db):
    gout_id = fw_spec.get("gaussian_output_id", {}).get(key)
    run = run_db.runs. \
        find_one({"_id": ObjectId(gout_id)})
    proceed_keys = fw_spec.get("proceed", {})
    for k, v in proceed_keys.items():
        if run["output"].get(k, run["output"]["output"].get(k)) != v:
            raise ValueError(
                f"The condition for {k} is not met, Terminating"
            )
    return run


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
            raise FileNotFoundError(
                "Please provide the database configurations")
    if isinstance(input_db, dict):
        db = GaussianCalcDb(**input_db)
    else:
        db = GaussianCalcDb.from_db_file(input_db)

    return db


def recursive_signature_remove(d):
    if isinstance(d, dict):
        return {i: recursive_signature_remove(j)
                for i, j in d.items() if not i.startswith('@')}
    else:
        return d