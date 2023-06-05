# coding: utf-8


# Defines the Gaussian database class.

import ssl
import sys
import json
import logging
import datetime

from abc import abstractmethod

import pandas as pd

from pymongo import ASCENDING, MongoClient

from monty.serialization import loadfn

from pymatgen.core.structure import Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher

from mispr import __version__ as infrastructure_version
from mispr.gaussian.utilities.metadata import get_chem_schema

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger()
ch = logging.StreamHandler(stream=sys.stdout)
logger.addHandler(ch)
logger.setLevel(20)


class GaussianCalcDb:
    """
    Class to help manage database insertions of molecules and Gaussian
    calculations.
    """
    def __init__(
        self,
        host,
        port=None,
        name=None,
        username=None,
        password=None,
        uri_mode=False,
        **kwargs
    ):
        self.host = host
        self.db_name = name
        self.user = username
        self.password = password
        self.port = int(port) if port else None
        try:
            if uri_mode:
                print("uri_mode")
                self.connection = MongoClient(host)
                dbname = host.split("/")[-1].split("?")[
                    0
                ]  # parse URI to extract dbname
                self.db = self.connection[dbname]
                print(self.db.list_database_names())
            else:
                print("not uri_mode")
                self.connection = MongoClient(
                    self.host,
                    self.port,
                    ssl=kwargs.get("ssl", False),
                    ssl_ca_certs=kwargs.get("ssl_ca_certs"),
                    ssl_certfile=kwargs.get("ssl_certfile"),
                    ssl_keyfile=kwargs.get("ssl_keyfile"),
                    ssl_pem_passphrase=kwargs.get("ssl_pem_passphrase"),
                    ssl_cert_reqs=kwargs.get("ssl_cert_reqs", ssl.CERT_NONE),
                    username=username,
                    password=password,
                    authsource=kwargs.get("authsource"),
                )
                self.db = self.connection[self.db_name]
                print(self.db.list_database_names())
        except:
            logger.error("Mongodb connection failed")
            raise Exception
        try:
            if self.user:
                self.db.authenticate(
                    self.user, self.password, source=kwargs.get("authsource", None)
                )
        except:
            logger.error("Mongodb authentication failed")
            raise ValueError
        self.molecules = self.db["molecules"]
        self.functional_groups = self.db["functional_groups"]
        self.derived_molecules = self.db["derived_molecules"]
        self.runs = self.db["runs"]

        self.build_indexes()

    @abstractmethod
    def build_indexes(self, background=True):
        """
        Builds indexes for the database.

        Args:
            background: if ``True``, this index should be created in the
            background.
        """
        self.molecules.create_index("inchi", unique=False, background=background)
        self.molecules.create_index("smiles", unique=True, background=background)
        self.molecules.create_index(
            "formula_alphabetical", unique=False, background=background
        )
        self.molecules.create_index("chemsys", unique=False, background=background)
        self.functional_groups.create_index("name", unique=True, background=background)
        self.derived_molecules.create_index(
            [
                ("inchi", ASCENDING),
                ("smiles", ASCENDING),
                ("formula_alphabetical", ASCENDING),
                ("chemsys", ASCENDING),
            ],
            unique=False,
            background=background,
        )
        self.runs.create_index(
            [
                ("inchi", ASCENDING),
                ("smiles", ASCENDING),
                ("formula_alphabetical", ASCENDING),
                ("chemsys", ASCENDING),
                ("type", ASCENDING),
                ("functional", ASCENDING),
                ("basis", ASCENDING),
                ("phase", ASCENDING),
                ("tag", ASCENDING),
            ],
            unique=False,
            background=background,
        )

    def query_molecules(self, query):
        """
        Query the molecules collection.

        Args:
            query (dict): a query document that selects which documents
                to include in the result set; e.g. keys can be inchi,
                smiles, chemsys, etc.

        Returns:
            pd.DataFrame: a dataframe of documents that match the query
        """
        projection = {
            "inchi": 1,
            "smiles": 1,
            "formula": 1,
            "formula_pretty": 1,
            "formula_anonymous": 1,
            "formula_alphabetical": 1,
            "chemsys": 1,
            "nsites": 1,
            "nelements": 1,
            "is_ordered": 1,
            "is_valid": 1,
        }
        return pd.DataFrame(list(self.molecules.find(query, projection)))

    def insert_molecule(self, mol, update_duplicates=False):
        """
        Insert a molecule into the molecules collection.

        Args:
            mol (Molecule): a pymatgen Molecule object to insert
            update_duplicates (bool): if ``True``, update the existing
                molecule in the db with the new one
        """
        mol_dict = get_chem_schema(mol)
        # Check if mol is already in db
        result = self.molecules.find_one({"inchi": mol_dict["inchi"]})
        if result:
            logger.info("{} already in database".format(mol_dict["inchi"]))
        # If mol is not in db, checks if the same molecule exists with a
        # different representation
        if result is None:
            m = MoleculeMatcher()
            result_list = list(
                self.molecules.find(
                    {"formula_alphabetical": mol_dict["formula_alphabetical"]}
                )
            )
            if result_list:
                logger.info(
                    "{} already in {} documents".format(
                        mol_dict["inchi"], len(result_list)
                    )
                )
            for i in result_list:
                saved_mol = Molecule.from_dict(i)
                if m.fit(saved_mol, mol):
                    result = i
                    mol_dict["inchi"] = i["inchi"]
                    break
        # If update_duplicates is set to True, updates existing document with
        # new geometry keeping the old inchi representation
        if result and update_duplicates:
            logger.info("Updating duplicate {}".format(mol_dict["inchi"]))
        if result is None or update_duplicates:
            mol_dict["version"] = infrastructure_version
            mol_dict["last_updated"] = datetime.datetime.utcnow()
            self.molecules.update_one(
                {"inchi": mol_dict["inchi"]}, {"$set": mol_dict}, upsert=True
            )
            return mol_dict["inchi"]
        else:
            logger.info("Skipping duplicate {}".format(mol_dict["inchi"]))
            return mol_dict["inchi"]

    def retrieve_molecule(self, inchi):
        """
        Retrieve a molecule from the molecules collection.

        Args:
            inchi (str): the inchi representation of the molecule

        Returns:
            dict
        """
        return self.molecules.find_one({"inchi": inchi})

    def delete_molecule(self, inchi):
        """
        Delete a molecule from the molecules collection.

        Args:
            inchi (str): the inchi representation of the molecule to delete
        """
        return self.molecules.delete_one({"inchi": inchi})

    def insert_run(self, grun):
        """
        Insert a Gaussian run into the runs collection.

        Args:
            grun (dict): a dictionary containing the Gaussian run
        """
        if "_id" in grun:
            stored_run = self.retrieve_run(_id=grun["_id"])
            if stored_run:
                return grun["_id"]
        grun["last_updated"] = datetime.datetime.utcnow()
        result = self.runs.insert_one(grun, bypass_document_validation=True)
        return result.inserted_id

    def retrieve_doc(
        self,
        collection_name,
        inchi=None,
        smiles=None,
        functional=None,
        basis=None,
        **kwargs
    ):
        """
        Retrieve a document from any collection of the database.

        Args:
            collection_name (str): the name of the collection, e.g.
                bde, molecules, runs, etc.
            inchi (str): the inchi representation of the molecule
            smiles (str): the smiles representation of the molecule
            functional (str): the name of the density functional
            basis (str): the name of the basis set
            **kwargs: other kwargs that can be used to query the collection

        Returns:
            list: a list of documents that match the query
        """
        query = {}
        if inchi:
            query["inchi"] = inchi
        if smiles:
            query["smiles"] = smiles
        if functional:
            query["functional"] = functional
        if basis:
            query["basis"] = basis
        query = {**query, **kwargs}
        return list(self.db[collection_name].find(query))

    def retrieve_run(
        self, inchi=None, smiles=None, functional=None, basis=None, **kwargs
    ):
        """
        Retrieve a run from the runs collection.

        Args:
            inchi (str): the inchi representation of the molecule
            smiles (str): the smiles representation of the molecule
            functional (str): the name of the density functional
            basis (str): the name of the basis set
            **kwargs: other kwargs that can be used to query the collection

        Returns:
            list: a list of documents that match the query
        """
        result = self.retrieve_doc(
            "runs",
            inchi=inchi,
            smiles=smiles,
            functional=functional,
            basis=basis,
            **kwargs
        )
        return result

    def move_runs(
        self,
        new_collection,
        inchi=None,
        smiles=None,
        functional=None,
        basis=None,
        **kwargs
    ):
        """
        Move documents from the runs collection to another collection.

        Args:
            new_collection (str): the name of the collection to move the runs to
            inchi (str): the inchi representation of the molecule
            smiles (str): the smiles representation of the molecule
            functional (str): the name of the density functional
            basis (str): the name of the basis set
            **kwargs : other kwargs that can be used to query the collection
        """
        runs = self.retrieve_run(inchi, smiles, functional, basis, **kwargs)
        self.db[new_collection].insert_many(runs)

    def update_run(
        self,
        new_values,
        inchi=None,
        smiles=None,
        job_type=None,
        functional=None,
        basis=None,
        phase=None,
        **kwargs
    ):
        """
        Update a document in the runs collection. If multiple documents
        match the query criteria, will select the first one.

        Args:
            new_values (dict): the new output values to update the
                document with, e.g. {"polarizability": 3.5}
            inchi (str): the inchi representation of the molecule
            smiles (str): the smiles representation of the molecule
            job_type (str): the type of job, e.g. "opt", "freq"
            functional (str): the name of the density functional
            basis (str): the name of the basis set
            phase (str): the phase of the job, e.g. "gas", "solution"
            **kwargs: other kwargs that can be used to query the collection
        """
        query = {}
        if inchi:
            query["inchi"] = inchi
        if smiles:
            query["smiles"] = smiles
        if job_type:
            query["type"] = job_type
        if functional:
            query["functional"] = functional
        if basis:
            query["basis"] = basis
        if phase:
            query["phase"] = phase
        query = {**query, **kwargs}
        run_ = self.retrieve_run(inchi, smiles, functional, basis, **kwargs)[0]
        run_["output"]["output"].update(new_values)
        self.runs.update_one(query, {"$set": run_})

    @classmethod
    def from_db_file(cls, db_file, admin=True):
        """
        Create a new database object from a database file.

        Args:
            db_file (str): the path to the database file
            admin (bool): whether to use admin credentials

        Returns:
            GaussianCalcDb
        """
        creds = loadfn(db_file)
        
        kwargs = creds.get(
            "mongoclient_kwargs", {}
        )  # any other MongoClient kwargs can go here ...
        
        if creds.get("uri_mode", False):
            return cls(
                creds["host"],
                uri_mode=True,
                **kwargs
            )
        else:
            if admin and "admin_user" not in creds and "readonly_user" in creds:
                raise ValueError(
                    "Trying to use admin credentials, "
                    "but no admin credentials are defined. "
                    "Use admin=False if only read_only "
                    "credentials are available."
                )

            if admin:
                user = creds.get("admin_user")
                password = creds.get("admin_password")
            else:
                user = creds.get("readonly_user")
                password = creds.get("readonly_password")

            if "authsource" in creds:
                kwargs["authsource"] = creds["authsource"]
            else:
                kwargs["authsource"] = creds["database"]

            return cls(
                creds["host"],
                int(creds.get("port", 27017)),
                creds["database"],
                user,
                password,
                **kwargs
            )

    def insert_fg(self, fg_file):
        """
        Insert functional groups into their collection using a json
        file. The file can contain one or more functional groups.

        Args:
            fg_file (str): the path to the json file
        """
        # file can contain one fg dictionary or multiple ones
        with open(fg_file) as f:
            func_group = json.load(f)
        list_of_groups = []
        search_query = []
        for i, j in func_group.items():
            list_of_groups.append({**j, "name": i})
            search_query.append(i)
        existing = self.functional_groups.find(
            {"name": {"$in": search_query}}, {"name": 1}
        )
        existing = set([i["name"] for i in existing])
        list_of_groups = [i for i in list_of_groups if i["name"] not in existing]
        self.functional_groups.insert_many(list_of_groups)

    def retrieve_fg(self, name):
        """
        Retrieve a functional group from the functional_groups collection.

        Args:
            name (str): the name of the functional group to retrieve

        Returns:
            dict
        """
        return self.functional_groups.find_one({"name": name})

    def insert_derived_mol(self, derived_mol, update_duplicates):
        """
        Insert a derived molecule into the derived_molecules collection.

        Args:
            derived_mol (Molecule): pymatgen.Molecule object
            update_duplicates (bool): whether to update duplicates if
                the molecule already exists
        """
        if isinstance(derived_mol, Molecule):
            derived_mol_dict = get_chem_schema(derived_mol)
            result = self.derived_molecules.find_one(
                {"inchi": derived_mol_dict["inchi"]}
            )
            if result:
                logger.info("{} already in database".format(derived_mol_dict["inchi"]))
            if result and update_duplicates:
                logger.info("Updating duplicate {}".format(derived_mol_dict["inchi"]))
            if result is None or update_duplicates:
                derived_mol_dict["last_updated"] = datetime.datetime.utcnow()
                self.derived_molecules.update_one(
                    {"inchi": derived_mol_dict["inchi"]},
                    {"$set": derived_mol_dict},
                    upsert=True,
                )
                return derived_mol_dict["inchi"]
        else:
            logger.info("No molecule provided")

    def insert_property(self, collection_name, property_dict, index, **kwargs):
        """
        Insert a document into a property collection in the database.
        If the property already exists,

        Args:
            collection_name (str): the name of the collection to insert
                the property into; e.g. bde, binding_energy, etc.
            property_dict (dict): the property dictionary to insert
            index (str, list[tuple]): the indexes to use for fast lookup;
            **kwargs: additional kwargs to pass to
                pymongo.collection.create_index

        Returns:

        """
        collection = self.db[collection_name]
        collection.create_index(index, **kwargs)
        collection.insert_one(property_dict)
