# coding: utf-8


# Defines firetasks for parsing Gaussian output files.

import os
import sys
import json
import logging
import datetime
import subprocess

from configparser import ConfigParser

import numpy as np
import pymongo
import networkx as nx
import matplotlib.image as img
import matplotlib.pyplot as plt

from bson.objectid import ObjectId

from fireworks.fw_config import CONFIG_FILE_DIR
from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize
from fireworks.utilities.fw_serializers import DATETIME_HANDLER

from pymatgen.io.gaussian import GaussianInput
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from mispr import __version__ as mispr_version
from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.utilities.gout import process_run
from mispr.gaussian.utilities.misc import pass_gout_dict
from mispr.gaussian.utilities.dbdoc import add_solvent_to_prop_dict
from mispr.gaussian.utilities.files import bibtex_parser
from mispr.gaussian.utilities.rdkit import (
    get_rdkit_mol,
    draw_rdkit_mol_with_highlighted_bonds,
)
from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.gaussian.utilities.db_utilities import get_db

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

DEFAULT_KEY = "gout_key"
HARTREE_TO_EV = 27.2114
HARTREE_TO_KJ = 2600
FARAD = 96.5
R = 8.31446


@explicit_serialize
class ProcessRun(FiretaskBase):
    required_params = ["run"]
    optional_params = [
        "operation_type",
        "db",
        "save_to_db",
        "save_to_file",
        "filename",
        "input_file",
        "gout_key",
        "format_chk",
        "formchk_cmd",
    ]

    def run_task(self, fw_spec):
        run = self["run"]
        operation_type = self.get("operation_type", "get_from_gout")
        input_file = self.get("input_file")
        working_dir = os.getcwd()
        db = self.get("db")

        if self.get("format_chk"):
            found_chk = False
            for file in os.listdir(working_dir):
                if file.endswith(".chk"):
                    found_chk = True
                    chk_file = file.strip(".chk")
                    cmd = self.get("formchk_cmd")
                    if not cmd:
                        cfg = ConfigParser()
                        cfg.read(CONFIG_FILE_DIR + "/config.ini")
                        cmd = cfg["RunCalc"]["formchkcmd"]
                    cmd = cmd.replace("$input_path$", file).replace(
                        "$output_path$", f"{chk_file}.fchk"
                    )

                    logger.info("Running command: {}".format(cmd))
                    logger.info("Running command: {}".format(cmd))
                    return_code = subprocess.call(cmd, shell=True)
                    logger.info(
                        "Finished running with return code: {}".format(return_code)
                    )
                if found_chk:
                    break
            if not found_chk:
                logger.info(f"No checkpoint file found in {working_dir}")

        gout_dict = process_run(
            operation_type=operation_type,
            run=run,
            input_file=input_file,
            working_dir=working_dir,
            db=db,
        )

        # get initial and final molecule
        input_mol = Molecule.from_dict(gout_dict["input"]["molecule"])
        output_mol = Molecule.from_dict(gout_dict["output"]["output"]["molecule"])

        # Build initial and final molgraphs
        input_molgraph = MoleculeGraph.with_local_env_strategy(input_mol, OpenBabelNN())
        output_molgraph = MoleculeGraph.with_local_env_strategy(
            output_mol, OpenBabelNN()
        )

        # Detect any structural change that occurred during calculation
        if input_molgraph.isomorphic_to(output_molgraph):
            change = "no_change"
        else:
            input_graph = input_molgraph.graph
            output_graph = output_molgraph.graph
            if nx.is_connected(input_graph.to_undirected()) and not nx.is_connected(
                output_graph.to_undirected()
            ):
                change = "unconnected_fragments"
            elif output_graph.number_of_edges() < input_graph.number_of_edges():
                change = "fewer_bonds"
            elif output_graph.number_of_edges() > input_graph.number_of_edges():
                change = "more_bonds"
            else:
                change = "bond_change"
        gout_dict["structural_change"] = change

        if "_id" in gout_dict:
            gout_dict["_id"] = ObjectId(gout_dict["_id"])
        if "tag" in fw_spec:
            gout_dict["tag"] = fw_spec["tag"]
        if "run_time" in fw_spec:
            gout_dict["wall_time (s)"] = fw_spec["run_time"]
        gout_dict["version"] = mispr_version
        if gout_dict["output"]["has_gaussian_completed"]:
            run_list = {}
            if self.get("save_to_db"):
                runs_db = get_db(db)
                run_id = runs_db.insert_run(gout_dict)
                run_list["run_id_list"] = run_id
                logger.info("Saved parsed output to db")

            if self.get("save_to_file"):
                filename = self.get("filename", "run")
                file = os.path.join(working_dir, f"{filename}.json")
                with open(file, "w") as f:
                    f.write(json.dumps(gout_dict, default=DATETIME_HANDLER))
                run_list["run_loc_list"] = file
                logger.info("Saved parsed output to json file")

            uid = self.get("gout_key")
            set_dict = {f"gaussian_output->{DEFAULT_KEY}": gout_dict}
            if uid:
                set_dict[f"gaussian_output->{uid}"] = gout_dict
            # fw_spec = {'gaussian_output: {DEFAULT_KEY: gout_dict, uid: gout_dict}}
            mod_dict = {"_set": set_dict}
            if run_list:
                mod_dict.update({"_push": run_list})
            return FWAction(mod_spec=mod_dict, propagate=True)
        else:
            raise ValueError(f"Gaussian did not complete normally, Terminating")


@explicit_serialize
class RetrieveGaussianOutput(FiretaskBase):
    """
    Returns a Gaussian output object from the database and converts it to a
    Gaussian input object
    """

    required_params = []
    optional_params = [
        "db",
        "gaussian_input_params",
        "run_id",
        "inchi",
        "functional",
        "basis",
        "type",
        "phase",
        "tag",
    ]

    def run_task(self, fw_spec):

        # if a Gaussian output dict is passed through fw_spec
        if fw_spec.get("gaussian_output"):
            run = pass_gout_dict(fw_spec, DEFAULT_KEY)
            # run = fw_spec['gaussian_output'][DEFAULT_KEY]

        elif self.get("run_id"):
            run = process_run(
                operation_type="get_from_run_id",
                run=self.get("run_id"),
                db=self.get("db"),
            )

        # otherwise, try to retrieve gaussian output from db
        else:
            try:
                query = {
                    "inchi": self.get("inchi"),
                    "smiles": self.get("smiles"),
                    "type": self.get("type"),
                    "functional": self.get("functional"),
                    "basis": self.get("basis"),
                    "phase": self.get("phase"),
                }
                if "tag" in self:
                    query["tag"] = self["tag"]
                run = process_run(
                    operation_type="get_from_run_query", run=query, db=self.get("db")
                )
            except pymongo.errors.ConnectionFailure as e:
                print("Could not connect to server: %s" % e)
            except Exception as e:
                raise ValueError(e)

        # create a gaussian input object from run
        if self.get("gaussian_input_params") is None:
            logger.info(
                "No gaussian input parameters provided; will use " "run parameters"
            )
        inputs = {}
        for k, v in run["input"].items():
            # use gaussian_input_params if defined, otherwise use run parameters
            inputs[f"{k}"] = self.get("gaussian_input_params", {}).get(
                f"{k}", run["input"].get(f"{k}")
            )
        inputs["molecule"] = run["output"]["output"]["molecule"]
        # TODO: check if this works in all cases
        #  (if we are saving to db and removing the class)
        # inputs['molecule'] = process_mol('get_from_run_dict', run)
        gaussin = GaussianInput.from_dict(inputs)
        fw_spec["gaussian_input"] = gaussin


@explicit_serialize
class ESPtoDB(FiretaskBase):
    required_params = ["keys"]
    optional_params = [
        "db",
        "save_to_db",
        "save_to_file",
        "additional_prop_doc_fields",
        "solvent_gaussian_inputs",
        "solvent_properties",
    ]

    def run_task(self, fw_spec):
        keys = self["keys"]
        db = self.get("db")
        gout_dict = [pass_gout_dict(fw_spec, i) for i in keys]
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = gout_dict + [pass_gout_dict(fw_spec, i + "_opt") for i in keys]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
        molecule = process_mol(
            "get_from_run_dict", gout_dict[-1], charge=gout_dict[-1]["input"]["charge"]
        )

        mol_schema = get_chem_schema(molecule)
        phase = gout_dict[-1]["phase"]

        # if one calculation is skipped, wall time is considered zero
        run_time = sum([gout.get("wall_time (s)", 0) for gout in full_gout_dict])

        esp_dict = {
            "molecule": molecule.as_dict(),
            "smiles": mol_schema["smiles"],
            "inchi": mol_schema["inchi"],
            "formula_alphabetical": mol_schema["formula_alphabetical"],
            "chemsys": mol_schema["chemsys"],
            "energy": gout_dict[-1]["output"]["output"]["final_energy"],
            "esp": gout_dict[-1]["output"]["output"]["ESP_charges"],
            "dipole_moment": gout_dict[-1]["output"]["output"]["dipole_moment"],
            "functional": gout_dict[-1]["functional"],
            "basis": gout_dict[-1]["basis"],
            "phase": phase,
            "tag": gout_dict[-1]["tag"],
            "state": "successful",
            "wall_time (s)": run_time,
            "version": mispr_version,
            "gauss_version": gout_dict[-1]["gauss_version"],
            "last_updated": datetime.datetime.utcnow(),
        }

        if phase == "solution":
            solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
            solvent_properties = self.get("solvent_properties", {})
            esp_dict = add_solvent_to_prop_dict(
                esp_dict, solvent_gaussian_inputs, solvent_properties
            )

        # check if polarizability is available (from freq calc of esp workflow)
        if "polarizability" in gout_dict[-2]["output"]["output"]:
            esp_dict["polarizability"] = gout_dict[-2]["output"]["output"][
                "polarizability"
            ]

        if self.get("additional_prop_doc_fields"):
            esp_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            esp_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get("save_to_db"):
            db = get_db(db)
            db.insert_property(
                "esp",
                esp_dict,
                [
                    ("formula_alphabetical", 1),
                    ("smiles", 1),
                    ("inchi", 1),
                    ("chemsys", 1),
                    ("functional", 1),
                    ("basis", 1),
                    ("tag", 1),
                ],
            )

        if fw_spec.get("run_loc_list"):
            esp_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get("save_to_file"):
            if "run_ids" in esp_dict:
                del esp_dict["run_ids"]
            with open("esp.json", "w") as f:
                f.write(json.dumps(esp_dict, default=DATETIME_HANDLER))

        logger.info("esp calculation complete")
        return FWAction()


@explicit_serialize
class NMRtoDB(FiretaskBase):
    required_params = ["keys"]
    optional_params = [
        "db",
        "save_to_db",
        "save_to_file",
        "additional_prop_doc_fields",
        "solvent_gaussian_inputs",
        "solvent_properties",
    ]

    def run_task(self, fw_spec):
        keys = self["keys"]
        db = self.get("db")
        gout_dict = [pass_gout_dict(fw_spec, i) for i in keys]
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = gout_dict + [pass_gout_dict(fw_spec, i + "_opt") for i in keys]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
        molecule = process_mol(
            "get_from_run_dict", gout_dict[-1], charge=gout_dict[-1]["input"]["charge"]
        )

        mol_schema = get_chem_schema(molecule)
        phase = gout_dict[-1]["phase"]

        # if one calculation is skipped, wall time is considered zero
        run_time = sum([gout.get("wall_time (s)", 0) for gout in full_gout_dict])

        nmr_dict = {
            "molecule": molecule.as_dict(),
            "smiles": mol_schema["smiles"],
            "inchi": mol_schema["inchi"],
            "formula_alphabetical": mol_schema["formula_alphabetical"],
            "chemsys": mol_schema["chemsys"],
            "energy": gout_dict[-1]["output"]["output"]["final_energy"],
            "tensor": gout_dict[-1]["output"]["output"]["tensor"],
            "functional": gout_dict[-1]["functional"],
            "basis": gout_dict[-1]["basis"],
            "phase": phase,
            "tag": gout_dict[-1]["tag"],
            "state": "successful",
            "wall_time (s)": run_time,
            "version": mispr_version,
            "gauss_version": gout_dict[-1]["gauss_version"],
            "last_updated": datetime.datetime.utcnow(),
        }

        if phase == "solution":
            solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
            solvent_properties = self.get("solvent_properties", {})
            nmr_dict = add_solvent_to_prop_dict(
                nmr_dict, solvent_gaussian_inputs, solvent_properties
            )

        if self.get("additional_prop_doc_fields"):
            nmr_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            nmr_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get("save_to_db"):
            db = get_db(db)
            db.insert_property(
                "nmr",
                nmr_dict,
                [
                    ("formula_alphabetical", 1),
                    ("smiles", 1),
                    ("inchi", 1),
                    ("chemsys", 1),
                    ("functional", 1),
                    ("basis", 1),
                    ("tag", 1),
                ],
            )

        if fw_spec.get("run_loc_list"):
            nmr_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get("save_to_file"):
            if "run_ids" in nmr_dict:
                del nmr_dict["run_ids"]
            with open("nmr.json", "w") as f:
                f.write(json.dumps(nmr_dict, default=DATETIME_HANDLER))

        logger.info("nmr calculation complete")
        return FWAction()


@explicit_serialize
class BindingEnergytoDB(FiretaskBase):
    required_params = ["index", "keys"]
    optional_params = [
        "db",
        "save_to_db",
        "save_to_file",
        "additional_prop_doc_fields",
        "solvent_gaussian_inputs",
        "solvent_properties",
    ]

    def run_task(self, fw_spec):
        index = self["index"]
        keys = self["keys"]
        db = self.get("db")
        gout_dict = [pass_gout_dict(fw_spec, i) for i in keys]
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = gout_dict + [pass_gout_dict(fw_spec, i + "_opt") for i in keys]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
        molecules = [
            process_mol("get_from_run_dict", gout, charge=gout["input"]["charge"])
            for gout in gout_dict
        ]
        final_energies = [
            gout["output"]["output"]["final_energy"] for gout in gout_dict
        ]
        # be_key = 'binding_energy_{}_{}_eV'.format(
        #     str(molecules[0].species[index[0]]) + str(index[0]),
        #     str(molecules[1].species[index[1]]) + str(len(molecules[0]) +
        #                                               index[1]))

        be_value = (
            final_energies[2] - (final_energies[0] + final_energies[1])
        ) * HARTREE_TO_EV

        be = {
            "sites": (index[0], len(molecules[0]) + index[1]),
            "atoms": (
                str(molecules[0].species[index[0]]),
                str(molecules[1].species[index[1]]),
            ),
            "value": be_value,
        }
        mol_schema = get_chem_schema(molecules[2])
        phase = gout_dict[-1]["phase"]
        # if one calculation is skipped, wall time is considered zero
        run_time = sum([gout.get("wall_time (s)", 0) for gout in full_gout_dict])

        be_dict = {
            "molecule": molecules[2].as_dict(),
            "smiles": mol_schema["smiles"],
            "inchi": mol_schema["inchi"],
            "formula_alphabetical": mol_schema["formula_alphabetical"],
            "chemsys": mol_schema["chemsys"],
            "energy": final_energies[2],
            "be_eV": be,
            "functional": gout_dict[-1]["functional"],
            "basis": gout_dict[-1]["basis"],
            "phase": gout_dict[-1]["phase"],
            "tag": gout_dict[-1]["tag"],
            "state": "successful",
            "wall_time (s)": run_time,
            "version": mispr_version,
            "gauss_version": gout_dict[-1]["gauss_version"],
            "last_updated": datetime.datetime.utcnow(),
        }

        if phase == "solution":
            solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
            solvent_properties = self.get("solvent_properties", {})
            be_dict = add_solvent_to_prop_dict(
                be_dict, solvent_gaussian_inputs, solvent_properties
            )

        if self.get("additional_prop_doc_fields"):
            be_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            be_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get("save_to_db"):
            db = get_db(db)
            db.insert_property(
                "binding_energy",
                be_dict,
                [
                    ("formula_alphabetical", 1),
                    ("smiles", 1),
                    ("inchi", 1),
                    ("chemsys", 1),
                    ("functional", 1),
                    ("basis", 1),
                    ("tag", 1),
                ],
            )

        if fw_spec.get("run_loc_list"):
            be_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get("save_to_file"):
            if "run_ids" in be_dict:
                del be_dict["run_ids"]
            with open("binding_energy.json", "w") as f:
                f.write(json.dumps(be_dict, default=DATETIME_HANDLER))

        logger.info("binding energy calculation complete")
        return FWAction()


@explicit_serialize
class IPEAtoDB(FiretaskBase):
    required_params = [
        "num_electrons",
        "states",
        "phases",
        "steps",
        "root_node_key",
        "keys",
        "pcet",
        "vertical",
    ]
    optional_params = [
        "db",
        "save_to_db",
        "save_to_file",
        "solvent_gaussian_inputs",
        "solvent_properties",
        "electrode_potentials",
        "additional_prop_doc_fields",
        "gibbs_elec",
        "gibbs_h",
    ]

    def run_task(self, fw_spec):
        db = self.get("db")
        num_electrons = self["num_electrons"]
        states = self["states"]
        phases = self["phases"]
        steps = self["steps"]
        root_node_key = self["root_node_key"]
        pcet = self["pcet"]
        keys = self["keys"]
        vertical = self["vertical"]
        gibbs_elec = self.get("gibbs_elec")
        gibbs_h = self.get("gibbs_h")

        data_dir = os.path.join(os.path.dirname(__file__), "../data")
        ref_potentials = {
            "hydrogen": {
                "potential": 4.43,
                "ref": bibtex_parser("h_pot.bib", data_dir),
            },
            "magnesium": {
                "potential": 2.375,
                "ref": bibtex_parser("mg_pot.bib", data_dir),
            },
            "lithium": {
                "potential": 1.40,
                "ref": bibtex_parser("li_pot.bib", data_dir),
            },
        }
        electrode_potentials = self.get("electrode_potentials") or {}
        electrode_potentials = {**ref_potentials, **electrode_potentials}
        gout_dict = {i: pass_gout_dict(fw_spec, i) for i in keys}
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = list(gout_dict.values()) + [
            pass_gout_dict(fw_spec, i + "_opt") for i in keys
        ]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
        molecule = process_mol(
            "get_from_run_dict",
            gout_dict[root_node_key],
            charge=gout_dict[root_node_key]["input"]["charge"],
        )
        final_energies = {
            i: j["output"]["output"]["final_energy"] for i, j in gout_dict.items()
        }
        free_energy_corrections = {
            i: j["output"]["output"]["corrections"]["Gibbs Free Energy"]
            for i, j in gout_dict.items()
        }
        free_energies = {
            key: final_energies[key] + free_energy_corrections[key]
            for key in gout_dict.keys()
        }

        # if one calculation is skipped, wall time is considered zero
        run_time = sum([gout.get("wall_time (s)", 0) for gout in full_gout_dict])

        def _name_(e1, e2, h1, h2):
            if pcet:
                return f"{e1}e{h1}h_{e2}e{h2}h"
            else:
                return f"{e1}e_{e2}e"

        def _gibbs_to_redox(s, g, e_count):
            redox_result = -s * g / ((1 / HARTREE_TO_KJ) * e_count * FARAD)
            return redox_result

        def _gibbs_to_pka(g):
            pka_result = (
                g * HARTREE_TO_KJ / (-1 / np.log10(np.exp(1)) * R / 1000 * 298.15)
            )
            return pka_result

        # TODO: PCET --> what to save (IP, EA, pKA, intermediate results)?
        # TODO: possible to have multiple H transfer in single step? same pka
        #  calc, divide by 2?
        # TODO: gibbs_h in gas vs sol or same?
        gibbs_dict = {}
        redox_dict = {}
        for state in states:
            process = "oxidation" if state == "cation" else "reduction"
            gibbs_dict[process] = {}
            redox_dict[process] = {}
            sign = 1 if state == "anion" else -1
            for phase in phases:
                gibbs_dict[process][phase] = {}
                redox_dict[process][phase] = {}
                # steps not involving h transfer
                if not pcet or sign > 0:
                    for h in range((num_electrons + 1) * pcet + (1 - pcet)):
                        n1 = f"{phase.lower()}_0_{h}"
                        n2 = f"{phase.lower()}_{sign * num_electrons}_{h}"
                        name = _name_(0, sign * num_electrons, h, h)
                        gibbs = (
                            free_energies[n2]
                            - free_energies[n1]
                            - sign * gibbs_elec * num_electrons
                        )
                        gibbs_dict[process][phase][name] = gibbs * HARTREE_TO_EV
                        redox = _gibbs_to_redox(sign, gibbs, num_electrons)
                        redox_dict[process][phase][name] = {}
                        redox_dict[process][phase][name]["raw"] = redox
                        for key, value in electrode_potentials.items():
                            redox_dict[process][phase][name][key] = (
                                redox - value["potential"]
                            )
                        if steps == "multi":
                            for i in range(num_electrons):
                                n1 = f"{phase.lower()}_{sign * i}_{h}"
                                n2 = f"{phase.lower()}_{sign * (i + 1)}_{h}"
                                name = _name_(sign * i, sign * (i + 1), h, h)
                                gibbs = (
                                    free_energies[n2]
                                    - free_energies[n1]
                                    - sign * gibbs_elec
                                )
                                gibbs_dict[process][phase][name] = gibbs * HARTREE_TO_EV
                                redox = _gibbs_to_redox(sign, gibbs, 1)
                                redox_dict[process][phase][name] = {}
                                redox_dict[process][phase][name]["raw"] = redox
                                for key, value in electrode_potentials.items():
                                    redox_dict[process][phase][name][key] = (
                                        redox - value["potential"]
                                    )
                # steps involving h transfer
                else:
                    # 1st ind --> elec, 2nd ind --> hydrogen
                    n1 = f"{phase.lower()}_0_0"
                    n2 = f"{phase.lower()}_{num_electrons}_{num_electrons}"
                    name = _name_(0, num_electrons, 0, num_electrons)
                    gibbs = (
                        free_energies[n2] - free_energies[n1] - gibbs_h * num_electrons
                    )
                    gibbs_dict[process][phase][name] = gibbs * HARTREE_TO_EV
                    for e in range(num_electrons + 1):
                        n1 = f"{phase.lower()}_{e}_0"
                        n2 = f"{phase.lower()}_{e}_{num_electrons}"
                        name = _name_(e, e, 0, num_electrons)
                        gibbs = (
                            free_energies[n2]
                            - free_energies[n1]
                            - gibbs_h * num_electrons
                        )
                        pka = _gibbs_to_pka(gibbs)
                        gibbs_dict[process][phase][name] = gibbs * HARTREE_TO_EV
                        redox_dict[process][phase][name] = pka
                        if steps == "multi":
                            for i in range(num_electrons):
                                n1 = f"{phase.lower()}_{e}_{i}"
                                n2 = f"{phase.lower()}_{e}_{i + 1}"
                                name = _name_(e, e, i, i + 1)
                                gibbs = free_energies[n2] - free_energies[n1] - gibbs_h
                                pka = _gibbs_to_pka(gibbs)
                                gibbs_dict[process][phase][name] = gibbs * HARTREE_TO_EV
                                redox_dict[process][phase][name] = pka

        mol_schema = get_chem_schema(molecule)
        ip_ea_dict = {
            "molecule": molecule.as_dict(),
            "smiles": mol_schema["smiles"],
            "inchi": mol_schema["inchi"],
            "formula_alphabetical": mol_schema["formula_alphabetical"],
            "chemsys": mol_schema["chemsys"],
            "num_electrons": num_electrons,
            "functional": gout_dict[root_node_key]["functional"],
            "basis": gout_dict[root_node_key]["basis"],
            "phase": phases,
            "steps": steps,
            "vertical": vertical,
            "pcet": pcet,
            "gibbs_elec": gibbs_elec * HARTREE_TO_EV,
            "gibbs_h": gibbs_h * HARTREE_TO_EV,
            "gibbs": gibbs_dict,
            "potentials": redox_dict,
            "standard_potentials": electrode_potentials,
            "tag": gout_dict[root_node_key]["tag"],
            "state": "successful",
            "wall_time (s)": run_time,
            "version": mispr_version,
            "gauss_version": gout_dict[root_node_key]["gauss_version"],
            "last_updated": datetime.datetime.utcnow(),
        }

        if "solution" in phases:
            solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
            solvent_properties = self.get("solvent_properties", {})
            ip_ea_dict = add_solvent_to_prop_dict(
                ip_ea_dict, solvent_gaussian_inputs, solvent_properties
            )

        if self.get("additional_prop_doc_fields"):
            ip_ea_dict.update(self.get("additional_prop_doc_fields"))

        if fw_spec.get("run_id_list"):
            ip_ea_dict["run_ids"] = fw_spec["run_id_list"]

        if self.get("save_to_db"):
            db = get_db(db)
            db.insert_property(
                "ip_ea",
                ip_ea_dict,
                [
                    ("formula_alphabetical", 1),
                    ("smiles", 1),
                    ("inchi", 1),
                    ("chemsys", 1),
                    ("functional", 1),
                    ("basis", 1),
                    ("tag", 1),
                ],
            )

        if fw_spec.get("run_loc_list"):
            ip_ea_dict["run_locs"] = fw_spec["run_loc_list"]
        if self.get("save_to_file"):
            if "run_ids" in ip_ea_dict:
                del ip_ea_dict["run_ids"]
            with open("ip_ea.json", "w") as f:
                f.write(json.dumps(ip_ea_dict, default=DATETIME_HANDLER))

        logger.info("ip/ea calculation complete")
        return FWAction()


@explicit_serialize
class BDEtoDB(FiretaskBase):
    required_params = []
    optional_params = [
        "principle_mol_key",
        "db",
        "save_to_db",
        "save_to_file",
        "solvent_gaussian_inputs",
        "solvent_properties",
        "additional_prop_doc_fields",
        "visualize",
        "color",
    ]

    def run_task(self, fw_spec):
        db = self.get("db")
        principle_mol_key = self.get("principle_mol_key", "ref_mol")
        keys = [principle_mol_key] + fw_spec["frag_keys"]
        bonds = fw_spec["bonds"]
        molecule_indices = fw_spec["molecule_indices"]
        gout_dict = {i: pass_gout_dict(fw_spec, i) for i in keys}
        # include opt fireworks in the list of gout_dict to calculate run time
        full_gout_dict = list(gout_dict.values()) + [
            pass_gout_dict(fw_spec, i + "_opt") for i in keys
        ]
        full_gout_dict = [i for i in full_gout_dict if i is not None]
        molecule = process_mol(
            "get_from_run_dict",
            gout_dict[principle_mol_key],
            charge=gout_dict[principle_mol_key]["input"]["charge"],
        )
        phase = gout_dict[principle_mol_key]["phase"]
        final_energies = {
            i: j["output"]["output"]["final_energy"] for i, j in gout_dict.items()
        }
        enthalpy_corrections = {
            i: j["output"]["output"]["corrections"]["Enthalpy"]
            for i, j in gout_dict.items()
        }
        enthalpies = {
            key: final_energies[key] + enthalpy_corrections[key]
            for key in gout_dict.keys()
        }

        run_time = sum([gout.get("wall_time (s)", 0) for gout in full_gout_dict])

        fragments = {}
        for gout_key, gout in gout_dict.items():
            if gout_key != principle_mol_key:
                frag = process_mol(
                    "get_from_run_dict", gout, charge=gout["input"]["charge"]
                )
                frag_schema = get_chem_schema(frag)
                fragments.update(
                    {
                        gout_key: {
                            "fragment": frag.as_dict(),
                            "smiles": frag_schema["smiles"],
                            "inchi": frag_schema["inchi"],
                            "formula_alphabetical": frag_schema["formula_alphabetical"],
                            "chemsys": frag_schema["chemsys"],
                            "energy": final_energies[gout_key],
                        }
                    }
                )

        bde_results = {}
        for ind, [bond, mol_ind] in enumerate(zip(bonds, molecule_indices)):
            frag_pairs = {}
            for i, j in enumerate(mol_ind):
                frag_ind = [j[count] for count in range(len(j))]
                frag_tuple = tuple([f"frag_{ind}" for ind in frag_ind])
                # try to calculate BDE from fragment pair results;
                # if the fragment is not found because the calculation failed,
                # move to the next pair and throw a warning
                try:
                    fragment_energies_sum = sum(
                        [enthalpies[frag] for frag in frag_tuple]
                    )
                    frag_pairs.update(
                        {
                            str(frag_ind): (
                                fragment_energies_sum - enthalpies[principle_mol_key]
                            )
                            * HARTREE_TO_EV
                        }
                    )
                except KeyError:
                    logger.warning("Results not found for pairs: {}".format(frag_tuple))
            bde_results.update({str(bond): frag_pairs})

        # if at least one BDE is calculated, continue with creating final dict
        if any(bde_results.values()):
            mol_schema = get_chem_schema(molecule)
            bde_dict = {
                "molecule": molecule.as_dict(),
                "smiles": mol_schema["smiles"],
                "inchi": mol_schema["inchi"],
                "formula_alphabetical": mol_schema["formula_alphabetical"],
                "chemsys": mol_schema["chemsys"],
                "energy": final_energies[principle_mol_key],
                "functional": gout_dict[principle_mol_key]["functional"],
                "basis": gout_dict[principle_mol_key]["basis"],
                "phase": phase,
                "fragments": fragments,
                "bde_eV": bde_results,
                "tag": gout_dict[principle_mol_key]["tag"],
                "state": "successful",
                "wall_time (s)": run_time,
                "version": mispr_version,
                "gauss_version": gout_dict[principle_mol_key]["gauss_version"],
                "last_updated": datetime.datetime.utcnow(),
            }

            if phase == "solution":
                solvent_gaussian_inputs = self.get("solvent_gaussian_inputs")
                solvent_properties = self.get("solvent_properties", {})
                bde_dict = add_solvent_to_prop_dict(
                    bde_dict, solvent_gaussian_inputs, solvent_properties
                )

            if self.get("visualize"):
                # if RDKit is not installed, throw a warning message and proceed
                # normally; visualization will not be done
                try:
                    num_bonds = len(bonds)
                    rdkit_mol = get_rdkit_mol(molecule)
                    color = tuple(
                        self.get("color", (197 / 255, 237 / 255, 223 / 255, 1))
                    )
                    max_bond = 10  # max # of bonds above which to change sizes
                    fig_size = tuple(
                        [int(num_bonds / max_bond) * 4 + i for i in (8, 6)]
                    )
                    _ = draw_rdkit_mol_with_highlighted_bonds(
                        rdkit_mol,
                        bonds,
                        colors=[color] * num_bonds,
                        filename="mol_bonds.png",
                    )
                    fig, (ax, picture) = plt.subplots(
                        ncols=2,
                        figsize=fig_size,
                        gridspec_kw={"width_ratios": [2, 1.2]},
                    )

                    mol_image = img.imread("mol_bonds.png")
                    picture.imshow(mol_image)
                    picture.axis("off")
                    os.remove("mol_bonds.png")

                    margin = 0.2
                    width = 1 - 2 * margin / num_bonds
                    gap = num_bonds + margin + width
                    x_ticks = []
                    xtick_labels = []
                    for ind, bond in enumerate(bonds):
                        data = [
                            gap + (1 + i) * width
                            for i in range(len(bde_results[str(bond)]))
                        ]
                        gap = data[-1] + width + int(num_bonds / max_bond)
                        x_ticks.append(sum(data) / len(data))
                        # xtick_labels.append(
                        #     '{}:{}-{}'.format(str(ind),
                        #                       (str(molecule.species[bond[0]])),
                        #                       (str(molecule.species[bond[1]]))))
                        xtick_labels.append("{}".format(str(ind)))
                        ax.bar(
                            data,
                            list(bde_results[str(bond)].values()),
                            width,
                            color=color,
                            edgecolor="black",
                            align="center",
                        )
                    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable="box")
                    ax.set_xticks(np.array(x_ticks))
                    # ax.set_xticklabels(xtick_labels, rotation=45)
                    ax.set_xticklabels(xtick_labels)
                    ax.tick_params(axis="x", length=0, labelsize=16)
                    ax.tick_params(axis="y", direction="in", length=8, labelsize=16)
                    ax.set_xlabel("Bond number", fontsize=16)
                    ax.set_ylabel("BDE (eV)", fontsize=16)
                    plt.savefig("bde.png", bbox_inches="tight", pad_inches=0.2)
                except Exception as e:
                    logger.warning(e)

            if self.get("additional_prop_doc_fields"):
                bde_dict.update(self.get("additional_prop_doc_fields"))

            if fw_spec.get("run_id_list"):
                bde_dict["run_ids"] = fw_spec["run_id_list"]

            if self.get("save_to_db"):
                db = get_db(db)
                db.insert_property(
                    "bde",
                    bde_dict,
                    [
                        ("formula_alphabetical", 1),
                        ("smiles", 1),
                        ("inchi", 1),
                        ("chemsys", 1),
                        ("functional", 1),
                        ("basis", 1),
                        ("tag", 1),
                    ],
                )

            if fw_spec.get("run_loc_list"):
                bde_dict["run_locs"] = fw_spec["run_loc_list"]
            if self.get("save_to_file"):
                if "run_ids" in bde_dict:
                    del bde_dict["run_ids"]
                with open("bde.json", "w") as f:
                    f.write(json.dumps(bde_dict, default=DATETIME_HANDLER))

            logger.info("bde calculation complete")

        else:
            logger.error("No BDE was calculated. Exiting ...")
            sys.exit()
        return FWAction()
