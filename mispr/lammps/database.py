# coding: utf-8

# Defines classes for handling mongodb interactions.

import sys
import logging

import datetime

from abc import abstractmethod

import pandas as pd

from pymongo import ASCENDING, MongoClient

from pymatgen.core.structure import Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher

from mispr.lammps.utilities.utils import process_ff_doc

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
            [("species_smiles", ASCENDING), ("force_field_types", ASCENDING)],
            unique=False,
            background=background,
        )
        self.runs.create_index(
            [
                ("species_smiles", ASCENDING),
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
        mol = parameter_dict["Molecule"]
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


if __name__ == "__main__":
    import os
    import collections as coll

    db_uri = "mongodb+srv://mbliss01:idlewide@gettingstarted.dt0sv.mongodb.net/lammps"

    test_connect = LammpsSysDb(db_uri, uri_mode=True)

    file_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "tests", "test_files", "data"
    )
    wat_mol = Molecule.from_file(os.path.join(file_dir, "SPC_E.pdb"))
    wat_param = {
        "Molecule": wat_mol,
        "Labels": ["ow", "hw", "hw"],
        "Masses": coll.OrderedDict({"ow": 16.000, "hw": 1.008}),
        "Nonbond": [[0.155394259, 3.16555789], [0.0, 0.0]],
        "Bonds": [{"coeffs": [553.0, 1], "types": [("ow", "hw")]}],
        "Angles": [{"coeffs": [100.0, 109.47], "types": [("hw", "ow", "hw")]}],
        "Dihedrals": [],
        "Impropers": [],
        "Improper Topologies": None,
        "Charges": [-0.8476, 0.4238, 0.4238],
        "common_name": "spc/e",
        "method": "literature",
        "doi": "https://doi.org/10.1021/j100308a038",
    }
    test_connect.insert_force_field(
        wat_param, wat_param["method"], doi=wat_param["doi"]
    )

    na_mol = Molecule.from_file(os.path.join(file_dir, "Na.pdb"))
    na_mol.set_charge_and_spin(1, 1)
    na_param = {
        "Molecule": na_mol,
        "Labels": ["na"],
        "Masses": coll.OrderedDict({"na": 22.99}),
        "Nonbond": [[0.02639002, 2.590733472]],
        "Bonds": [],
        "Angles": [],
        "Dihedrals": [],
        "Impropers": [],
        "Improper Topologies": None,
        "Charges": [1.0],
        "method": "literature",
        "method_note": "spc/e",
        "doi": "https://doi.org/10.1021/ct500918t",
    }

    test_connect.insert_force_field(na_param, na_param["method"], na_param["doi"], True)
