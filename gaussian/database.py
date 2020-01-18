import datetime
from abc import ABCMeta, abstractmethod
from logging import getLogger, INFO
from pymatgen.io.babel import BabelMolAdaptor
import pybel as pb
from pymongo import MongoClient, ReturnDocument, ASCENDING


logger = getLogger(__name__)
logger.setLevel(INFO)


class GaussianCalcDb:

    def __init__(self, host, port=None, name=None, username=None, password=None,
                 uri_mode=False, **kwargs):
        self.host = host
        self.db_name = name
        self.user = username
        self.password = password
        self.port = int(port) if port else None

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
        self.runs.create_index([("smiles", ASCENDING),
                                ("type", ASCENDING),
                                ("functional", ASCENDING),
                                ("basis", ASCENDING)],
                               unique=False, background=background)

    def insert_molecule(self, molecule, update_duplicates=True):
        """
        Insert the task document ot the database collection.
        Args:
            d (dict): task document
            update_duplicates (bool): whether to update the duplicates
        """
        mol_smiles = GaussianCalcDb.get_smiles(molecule)
        molecule_dict = molecule.as_dict()
        result = self.molecules.find_one({"smiles": mol_smiles})
        if result is None or update_duplicates:
            molecule_dict["last_updated"] = datetime.datetime.utcnow()
            self.molecules.update_one({"smiles": mol_smiles},
                                       {"$set": molecule_dict}, upsert=True)
            return mol_smiles
        else:
            logger.info("Skipping duplicate {}".format(mol_smiles))
            return None

    def insert_run(self, run, update_duplicates=True):
        """
        Insert the task document to the database collection.
        Args:
            d (dict): task document
            update_duplicates (bool): whether to update the duplicates
        """
        query = {i:j for i, j in run.items() if i not in ['input', 'output']}
        result = self.runs.find_one(query)
        if result is None or update_duplicates:
            run["last_updated"] = datetime.datetime.utcnow()
            self.runs.update_one(query, {"$set": run}, upsert=True)
            return query
        else:
            logger.info("Skipping duplicate")
            return None
