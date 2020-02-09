import datetime
from abc import ABCMeta, abstractmethod
import logging
import sys

from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.core.structure import Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher

import pybel as pb
from pymongo import MongoClient, ReturnDocument, ASCENDING

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
                                              ssl=kwargs.get('ssl'),
                                              ssl_ca_certs=kwargs.get('ssl_ca_certs'),
                                              ssl_certfile=kwargs.get('ssl_certfile'),
                                              ssl_keyfile=kwargs.get('ssl_keyfile'),
                                              ssl_pem_passphrase=kwargs.get('ssl_pem_passphrase'),
                                              username=kwargs.get('user'),
                                              password=kwargs.get('password'),
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

    @staticmethod
    def get_smiles(mol):
        """
        Returns canonical SMILES representation of a pymatgen molecule object.
        """
        a = BabelMolAdaptor(mol)
        pm = pb.Molecule(a.openbabel_mol)
        return pm.write("can").strip()

    @abstractmethod
    def build_indexes(self, indexes=None, background=True):
        """
         Build the indexes.
         Args:
             indexes (list): list of single field indexes to be built.
             background (bool): Run in the background or not.
         """
        self.molecules.create_index("smiles", unique=True, background=background)
        self.molecules.create_index("formula", unique=False,
                                    background=background)
        self.runs.create_index([("smiles", ASCENDING),
                                ("type", ASCENDING),
                                ("functional", ASCENDING),
                                ("basis", ASCENDING)],
                               unique=False, background=background)

    def insert_molecule(self, molecule, update_duplicates=False):
        """
        Insert the task document ot the database collection. Does not allow
        inserting an existing molecule in a new document.
        Args:
            d (dict): task document
            update_duplicates (bool): whether to update the duplicates
        """
        mol_smiles = GaussianCalcDb.get_smiles(molecule)
        molecule_dict = molecule.as_dict()
        molecule_dict['formula'] = molecule.composition.formula
        # Check if smile is already in db
        result = self.molecules.find_one({"smiles": mol_smiles})
        if result:
            logger.info("{} already in database".format(mol_smiles))
        # If smile is not in db, checks if the same molecule exists with a
        # different smile representation
        if result is None:
            m = MoleculeMatcher()
            result_list = \
                self.molecules.find({'formula': molecule_dict['formula']})
            if result_list:
                logger.info("{} already in {} documents".format(mol_smiles,
                                                                len(result_list)))
            for i in result_list:
                saved_mol = Molecule.from_dict(i)
                if m.fit(saved_mol, molecule):
                    result = i
                    mol_smiles = i['smiles']
                    break
        # If update_duplicates is set to True, updates existing document with
        # new geometry keeping the old smiles representation
        if result and update_duplicates:
            logger.info("Updating duplicate {}".format(mol_smiles))
        if result is None or update_duplicates:
            molecule_dict["last_updated"] = datetime.datetime.utcnow()
            self.molecules.update_one({"smiles": mol_smiles},
                                       {"$set": molecule_dict}, upsert=True)
            return mol_smiles
        else:
            logger.info("Skipping duplicate {}".format(mol_smiles))
            return mol_smiles

    def retrieve_molecule(self, smiles):
        return self.molecules.find_one({"smiles": smiles})

    def delete_molecule(self, smiles):
        return self.molecules.delete_one({"smiles": smiles})

    def insert_run(self, grun):
        grun["last_updated"] = datetime.datetime.utcnow()
        self.runs.insert_one(grun)
        # query = {i:j for i, j in grun.items() if i not in ['input', 'output']}
        # result = self.runs.find_one(query)
        # if result is None or update_duplicates:
        #     grun["last_updated"] = datetime.datetime.utcnow()
        #     self.runs.update_one(query, {"$set": grun}, upsert=True)
        #     return query
        # else:
        #     logger.info("Skipping duplicate")
        #     return None

    def retrieve_run(self, smiles=None, job_type=None, functional=None,
                     basis=None, **kwargs):
        query = {}
        if smiles:
            query['smiles'] = smiles
        if job_type:
            query['type'] = job_type
        if job_type:
            query['functional'] = functional
        if job_type:
            query['basis'] = basis
        query = {**query, **kwargs}
        print(query)
        return list(self.runs.find(query))

    def move_runs(self, new_collection, smiles=None, job_type=None,
                  functional=None, basis=None, **kwargs):
        runs = self.retrieve_run(smiles, job_type, functional, basis, **kwargs)
        self.db[new_collection].insert_many(runs)


    #TODO: check if a molecule has already been calculated in the database