# coding: utf-8

# Defines classes for handling mongodb interactions

from abc import abstractmethod
import logging
import sys

import pandas as pd
import json

from pymatgen.core.structure import Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher
from pymongo import MongoClient, ASCENDING

logger = logging.getLogger()
ch = logging.StreamHandler(stream=sys.stdout)
logger.addHandler(ch)
logger.setLevel(20)

class CustomForcefieldDb:

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
                dbname = host.split('/')[-1].split('?')[0] # parse URI to extract dbname
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
                                     source=kwargs.get('authsource', None))
        except:
            logger.error("Mongodb authentication failed")
            raise ValueError
        # self.force_fields = self.db['force_fields']
        self.systems = self.db['systems']
        self.runs = self.db['runs']

        self.build_indexes()

    @abstractmethod
    def build_indexes(self, indexes=None, background=True):
        self.systems.create_index([('species_smiles', ASCENDING),
                                   ('mixture_data', ASCENDING),
                                   ('box_size', ASCENDING),
                                   ('species_ff', ASCENDING)],
                                  unique=False, background=background)
        self.runs.create_index([('species_smiles', ASCENDING),
                                ('mixture_data', ASCENDING),
                                ('box_size', ASCENDING),
                                ('species_ff', ASCENDING),
                                ('job_type', ASCENDING)],
                               unique=False, background=background)

    def query_systems(self, query):
        projection = {'species_smiles': 1}
        pass