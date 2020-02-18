import datetime
from abc import abstractmethod
import logging
import sys
import pandas as pd

from pymatgen.core.structure import Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from monty.serialization import loadfn
from pymongo import MongoClient, ASCENDING
from infrastructure.gaussian.utils.utils import get_chem_schema

logger = logging.getLogger()
ch = logging.StreamHandler(stream=sys.stdout)
logger.addHandler(ch)
logger.setLevel(20)


class GaussianCalcDb:

    def __init__(self, host, port=None, name=None, username=None, password=None,
                 uri_mode=False, **kwargs):
        self.host = host
        self.db_name = name
        self.user = username
        self.password = password
        self.port = int(port) if port else None
        # TODO: read database credentials from a config file instead (atomate)
        try:
            if uri_mode:
                self.connection = MongoClient(host)
                dbname = host.split('/')[-1].split('?')[
                    0]  # parse URI to extract dbname
                self.db = self.connection[dbname]
            else:
                self.connection = MongoClient(self.host, self.port,
                                              ssl=kwargs.get('ssl', False),
                                              ssl_ca_certs=kwargs.get('ssl_ca_certs'),
                                              ssl_certfile=kwargs.get('ssl_certfile'),
                                              ssl_keyfile=kwargs.get('ssl_keyfile'),
                                              ssl_pem_passphrase=kwargs.get('ssl_pem_passphrase'),
                                              username=username,
                                              password=password,
                                              authsource=kwargs.get('authsource'))
                self.db = self.connection[self.db_name]
        except:
            logger.error("Mongodb connection failed")
            raise Exception
        try:
            if self.user:
                self.db.authenticate(self.user, self.password,
                                     source=kwargs.get("authsource", None))
        except:
            logger.error("Mongodb authentication failed")
            raise ValueError
        self.molecules = self.db['molecules']
        self.runs = self.db['runs']

        self.build_indexes()

    @abstractmethod
    def build_indexes(self, indexes=None, background=True):
        """
         Build the indexes.
         Args:
             indexes (list): list of single field indexes to be built.
             background (bool): Run in the background or not.
         """
        self.molecules.create_index("smiles", unique=True, background=background)
        self.molecules.create_index(
            "formula_alphabetical", unique=False, background=background)
        self.molecules.create_index(
            "chemsys", unique=False, background=background)
        self.runs.create_index([("smiles", ASCENDING),
                                ("type", ASCENDING),
                                ("functional", ASCENDING),
                                ("basis", ASCENDING)],
                               unique=False, background=background)

    def query_molecules(self, query):
        projection = {'smiles': 1,
                      "formula": 1,
                      "formula_pretty": 1,
                      "formula_anonymous": 1,
                      "formula_alphabetical": 1,
                      "chemsys": 1,
                      "nsites": 1,
                      "nelements": 1,
                      "is_ordered": 1,
                      "is_valid": 1}
        return pd.DataFrame(list(self.molecules.find(query, projection)))

    def insert_molecule(self, mol, update_duplicates=False):
        """
        Insert the task document ot the database collection. Does not allow
        inserting an existing molecule in a new document.
        Args:
            d (dict): task document
            update_duplicates (bool): whether to update the duplicates
        """
        mol_dict = get_chem_schema(mol)
        # Check if smile is already in db
        result = self.molecules.find_one({"smiles": mol_dict["smiles"]})
        if result:
            logger.info("{} already in database".format(mol_dict["smiles"]))
        # If smile is not in db, checks if the same molecule exists with a
        # different smile representation
        if result is None:
            m = MoleculeMatcher()
            result_list = list(self.molecules.
                               find({'formula_alphabetical': mol_dict["formula_alphabetical"]}))
            if result_list:
                logger.info("{} already in {} documents".format(mol_dict["smiles"],
                                                                len(result_list)))
            for i in result_list:
                saved_mol = Molecule.from_dict(i)
                if m.fit(saved_mol, mol):
                    result = i
                    mol_dict["smiles"] = i['smiles']
                    break
        # If update_duplicates is set to True, updates existing document with
        # new geometry keeping the old smiles representation
        if result and update_duplicates:
            logger.info("Updating duplicate {}".format(mol_dict["smiles"]))
        if result is None or update_duplicates:
            mol_dict["last_updated"] = datetime.datetime.utcnow()
            self.molecules.update_one({"smiles": mol_dict["smiles"]},
                                       {"$set": mol_dict}, upsert=True)
            return mol_dict["smiles"]
        else:
            logger.info("Skipping duplicate {}".format(mol_dict["smiles"]))
            return mol_dict["smiles"]

    def retrieve_molecule(self, smiles):
        return self.molecules.find_one({"smiles": smiles})

    def delete_molecule(self, smiles):
        return self.molecules.delete_one({"smiles": smiles})

    def insert_run(self, grun):
        # TODO: perform checks if a similar calculation is already in the db
        grun["last_updated"] = datetime.datetime.utcnow()
        result = self.runs.insert_one(grun, bypass_document_validation=True)
        return result.inserted_id

    def retrieve_run(self, smiles, job_type=None, functional=None,
                     basis=None, **kwargs):
        query = {}
        if smiles:
            query['smiles'] = smiles
        if job_type:
            query['type'] = job_type
        if functional:
            query['functional'] = functional
        if basis:
            query['basis'] = basis
        query = {**query, **kwargs}
        return list(self.runs.find(query))

    def move_runs(self, new_collection, smiles=None, job_type=None,
                  functional=None, basis=None, **kwargs):
        runs = self.retrieve_run(smiles, job_type, functional, basis, **kwargs)
        self.db[new_collection].insert_many(runs)


    #TODO: check if a molecule has already been calculated in the database