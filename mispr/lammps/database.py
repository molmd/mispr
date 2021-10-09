# coding: utf-8

# Defines classes for handling mongodb interactions.

import sys
import logging
import datetime

from abc import abstractmethod

import pandas as pd
from monty.serialization import loadfn

from pymongo import ASCENDING, MongoClient

from pymatgen.core.structure import Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher

from mispr.lammps.utilities.utilities import process_ff_doc

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"

logger = logging.getLogger()
ch = logging.StreamHandler(stream=sys.stdout)
logger.addHandler(ch)
logger.setLevel(20)


class LammpsSysDb:
    def __init__(
        self,
        host,
        port=None,
        name=None,
        username=None,
        password=None,
        uri_mode=False,
        **kwargs,
    ):
        self.host = host
        self.db_name = name
        self.user = username
        self.password = password
        self.port = int(port) if port else None
        try:
            if uri_mode:
                self.connection = MongoClient(host)
                # parse URI to extract dbname
                dbname = host.split("/")[-1].split("?")[0]
                self.db = self.connection[dbname]
            else:
                self.connection = MongoClient(
                    self.host,
                    self.port,
                    ssl=kwargs.get("ssl", False),
                    ssl_ca_certs=kwargs.get("ssl_ca_certs"),
                    ssl_certfile=kwargs.get("ssl_certfile"),
                    ssl_keyfile=kwargs.get("ssl_keyfile"),
                    ssl_pem_passphrase=kwargs.get("ssl_pem_passphrase"),
                    username=username,
                    password=password,
                    authsource=kwargs.get("authsource"),
                )
                self.db = self.connection[self.db_name]
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
        self.force_fields = self.db["force_fields"]
        self.systems = self.db["systems"]
        self.runs = self.db["runs"]

        self.build_indexes()

    @abstractmethod
    def build_indexes(self, indexes=None, background=True):
        self.force_fields.create_index(
            [("smiles", ASCENDING), ("method", ASCENDING), ("doi", ASCENDING)],
            unique=True,
            background=background,
        )
        self.systems.create_index(
            [("smiles", ASCENDING), ("force_field_types", ASCENDING)],
            unique=False,
            background=background,
        )
        self.runs.create_index(
            [
                ("smiles", ASCENDING),
                ("mixture_data", ASCENDING),
                ("box", ASCENDING),
                ("job_type", ASCENDING),
            ],
            unique=False,
            background=background,
        )

    def query_force_fields(self, query):
        projection = {
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
            "method": 1,
            "doi": 1,
        }
        return pd.DataFrame(list(self.force_fields.find(query, projection)))

    def insert_force_field(
        self, parameter_dict, method, doi=None, update_duplicates=False, **kwargs
    ):
        ff_doc = process_ff_doc(parameter_dict, method, doi, **kwargs)
        # mol = parameter_dict["Molecule"]
        mol = Molecule.from_dict(ff_doc)
        # Check if smiles is already in db
        result = self.force_fields.find_one(
            {"smiles": ff_doc["smiles"], "method": method, "doi": doi}
        )
        if result:
            logger.info(
                f"Parameters for {ff_doc['smiles']} using {method} "
                f"according to {doi} are already in database"
            )
        # If smile is not in db, checks if the same molecule exists with a
        # different smile representation
        if result is None:
            m = MoleculeMatcher()
            result_list = list(
                self.force_fields.find(
                    {
                        "formula_alphabetical": ff_doc["formula_alphabetical"],
                        "method": method,
                        "doi": doi,
                    }
                )
            )
            if result_list:
                logger.info(
                    f"Parameters for {ff_doc['smiles']} using "
                    f"{method} according to {doi} are already in "
                    f"{len(result_list)} documents"
                )
            for i in result_list:
                saved_mol = Molecule.from_dict(i)
                if m.fit(saved_mol, mol):
                    result = i
                    ff_doc["smiles"] = i["smiles"]
                    break
        # If update_duplicates is set to True, updates existing document with
        # new geometry keeping the old smiles representation
        if result and update_duplicates:
            logger.info(
                f"Updating duplicate parameters for "
                f"{ff_doc['smiles']} using {method} according to "
                f"{doi}"
            )
        if result is None or update_duplicates:
            ff_doc.update(parameter_dict)
            ff_doc["method"] = method
            ff_doc["doi"] = doi
            ff_doc["last_updated"] = datetime.datetime.utcnow()
            self.force_fields.update_one(
                {"smiles": ff_doc["smiles"], "method": method, "doi": doi},
                {"$set": ff_doc},
                upsert=True,
            )
            return ff_doc["smiles"]
        else:
            logger.info(
                f"Skipping duplicate parameters for "
                f"{ff_doc['smiles']} using {method} according to "
                f"{doi}"
            )
            return ff_doc["smiles"], method, doi

    def retrieve_force_field(self, smiles, method, doi=None):
        return self.force_fields.find_one(
            {"smiles": smiles, "method": method, "doi": doi}
        )

    def delete_force_field(self, smiles, method, doi=None):
        return self.force_fields.delete_one(
            {"smiles": smiles, "method": method, "doi": doi}
        )

    def insert_run(self, lmp_run):
        # TODO: perform checks if a similar calculation is already in the db
        if "_id" in lmp_run:
            stored_run = self.retrieve_run(_id=lmp_run["_id"])
            if stored_run:
                return lmp_run["_id"]
        lmp_run["last_updated"] = datetime.datetime.utcnow()
        result = self.runs.insert_one(lmp_run, bypass_document_validation=True)
        return result.inserted_id

    def retrieve_run(self, _id):
        return self.runs.find_one({"_id": _id})

    def insert_system(self, sys_doc):
        sys_doc["last_updated"] = datetime.datetime.utcnow()
        result = self.systems.insert_one(sys_doc, bypass_document_validation=True)
        return result.inserted_id

    @classmethod
    def from_db_file(cls, db_file, admin=True):
        # TODO: move this to a common database module (common with Gaussian)
        creds = loadfn(db_file)

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

        kwargs = creds.get(
            "mongoclient_kwargs", {}
        )  # any other MongoClient kwargs can go here ...

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
