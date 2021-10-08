# coding: utf-8


# Defines firetasks for reading output of Ambertools programs and LAMMPS.

import os
import json
import inspect

import numpy as np
import pandas as pd

from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from pymatgen.io.ambertools import PrmtopParser
from pymatgen.core.structure import Molecule

from mdproptools.structural.rdf_cn import calc_atomic_rdf, calc_molecular_rdf
from mdproptools.dynamical.diffusion import Diffusion

import mispr.lammps.utilities.utilities as iluu

from mispr.lammps.defaults import MSD_SETTINGS, RDF_SETTINGS, DIFF_SETTINGS
from mispr.gaussian.utilities.metadata import get_mol_formula

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"

GAFF_DOI = "https://doi.org/10.1002/jcc.20035"


@explicit_serialize
class ProcessPrmtop(FiretaskBase):
    _fw_name = "Process Prmtop"
    required_params = ["molecule"]
    optional_params = [
        "working_dir",
        "db",
        "prmtop_path",
        "prmtop_filename",
        "prmtop_dir",
        "unique_molecule_name",
        "system_force_field_dict",
        "doi",
        "save_ff_to_db",
        "save_ff_to_file",
        "ff_filename",
    ]

    def run_task(self, fw_spec):

        working_dir = self.get("working_dir", os.getcwd())
        os.chdir(working_dir)

        db = self.get("db", None)
        ff_filename = self.get("ff_filename", "ff.json")
        save_to_db = self.get("save_ff_to_db")
        save_to_file = self.get("save_ff_to_file")

        if isinstance(self.get("prmtop_path"), str):
            prmtop_file_path = self.get("prmtop_path")

        elif isinstance(self.get("prmtop_filename"), str):
            prmtop_filename = self.get("prmtop_filename")
            prmtop_dir = self.get("prmtop_dir", working_dir)
            prmtop_file_path = os.path.join(prmtop_dir, prmtop_filename)

        else:
            raise Exception(
                'Neither name nor path to a prmtop file was \
                            provided; either provide prmtop file name as \
                            "prmtop_filename" or path to prmtop as \
                            "prmtop_path".'
            )

        if not os.path.exists(prmtop_file_path):
            raise Exception(
                '"prmtop_file_path" is not a valid path; check \
                            that "prmtop_path" is correct or that \
                            "prmtop_filename" exists in "prmtop_dir"'
            )

        molecule = self.get("molecule")

        if not isinstance(molecule, Molecule):
            raise Exception('"molecule" is not a pymatgen Molecule object')

        # TODO: decide what unique id for the ff labels. If staying as smiles \
        #  string, decide how to get smiles.
        unique_mol_name = self.get("unique_molecule_name", get_mol_formula(molecule))
        ff_param_dict_general = PrmtopParser(prmtop_file_path, molecule, "").to_dict()
        ff_param_dict_system = iluu.add_ff_labels_to_dict(
            ff_param_dict_general, unique_mol_name
        )

        gaff_doi = self.get("gaff_doi", GAFF_DOI)

        if db and save_to_db:
            ff_db = iluu.get_db(input_db=db)
            ff_db.insert_force_field(ff_param_dict_general, "gaff", doi=gaff_doi)

        if save_to_file:
            ff_doc = iluu.process_ff_doc(ff_param_dict_general, "gaff", doi=gaff_doi)
            with open(os.path.join(working_dir, ff_filename), "w") as file:
                json.dump(ff_doc, file)

        sys_ff_dict = fw_spec.get(
            "system_force_field_dict", self.get("system_force_field_dict", {})
        )
        sys_ff_dict[unique_mol_name] = ff_param_dict_system

        return FWAction(
            stored_data={
                "ff_param_dict_system": ff_param_dict_system,
                "ff_param_dict_general": ff_param_dict_general,
            },
            update_spec={"system_force_field_dict": sys_ff_dict},
        )


@explicit_serialize
class CalcDiff(FiretaskBase):
    _fw_name = "Calculate Diffusion"
    required_params = []
    optional_params = ["working_dir", "diff_settings"]

    def run_task(self, fw_spec):
        working_dir = self.get("working_dir", os.getcwd())
        os.makedirs(working_dir, exist_ok=True)
        diff_settings = {**MSD_SETTINGS.copy(), **DIFF_SETTINGS.copy()}
        diff_settings.update(self.get("diff_settings", {}))
        msd_method = diff_settings["msd_method"]
        outputs_dir = diff_settings.get("outputs_dir",
                                    os.path.abspath(
                                        os.path.join(working_dir, "..", "nvt")))
        msd_df = None
        diff = Diffusion(
            timestep=diff_settings["timestep"],
            units=diff_settings["units"],
            outputs_dir=outputs_dir,
            diff_dir=working_dir,
        )
        if msd_method == "from_dump":
            num_mols = fw_spec.get("num_mols_list", None)
            num_mols = [int(i) for i in num_mols] if num_mols else None
            num_atoms_per_mol = fw_spec.get("num_atoms_per_mol", None)
            mass = fw_spec.get("default_masses", None)
            file_pattern = diff_settings.pop("file_pattern", "dump.nvt.*.dump")

            msd_df = diff.get_msd_from_dump(
                file_pattern,
                num_mols=num_mols,
                num_atoms_per_mol=num_atoms_per_mol,
                mass=mass,
                **{
                    i: j
                    for i, j in diff_settings.items()
                    if i in inspect.getfullargspec(diff.get_msd_from_dump).args
                }
            )
        elif msd_method == "from_log":
            file_pattern = diff_settings.get("file_pattern", "log.lammps*")
            msd_df = diff.get_msd_from_log(file_pattern)

        else:
            raise ValueError(f"Unsupported msd calculation method: {msd_method}"
                             f"Supported methods are from_dump and from_log")

        diffusion = diff.calc_diff(
            msd_df[0],
            **{
                i: j
                for i, j in diff_settings.items()
                if i in inspect.getfullargspec(diff.calc_diff).args
            }
        )

        if msd_method == "from_dump" and diff_settings["avg_interval"] and diff_settings["diff_dist"]:
            diff.get_diff_dist(msd_df[2], **{
                i: j
                for i, j in diff_settings.items()
                if i in inspect.getfullargspec(diff.get_diff_dist).args
            })

        smiles_list = fw_spec.get("smiles", [])
        n_mols_dict = fw_spec.get("nmols", {})
        num_mols_list = fw_spec.get("num_mols_list", [])
        lmp_box = fw_spec.get("box", None)
        num_atoms_per_mol_list = fw_spec.get("num_atoms_per_mol", [])
        default_masses_list = fw_spec.get("default_masses", [])
        recalc_masses_list = fw_spec.get("recalc_masses", [])

        return FWAction(
            update_spec={
                "diff_folder_path": working_dir,
                "smiles": smiles_list,
                "nmols": n_mols_dict,
                "num_mols_list": num_mols_list,
                "box": lmp_box,
                "num_atoms_per_mol": num_atoms_per_mol_list,
                "default_masses": default_masses_list,
                "recalc_masses": recalc_masses_list,
            }
        )


@explicit_serialize
class GetRDF(FiretaskBase):
    _fw_name = "Get RDF"
    required_params = []
    optional_params = ["rdf_settings", "working_dir"]

    def run_task(self, fw_spec):
        # Get inputs to rdf function from defaults or from user, which then
        # updates the defaults
        rdf_settings = RDF_SETTINGS.copy()
        rdf_settings.update(self.get("rdf_settings", {}))

        # atomic or molecular (for CoM)
        rdf_type = rdf_settings.get("rdf_type", "atomic")

        working_dir = self.get("working_dir", os.getcwd())
        os.makedirs(working_dir, exist_ok=True)

        use_default_atom_ids = rdf_settings.get("use_default_atom_ids", False)

        csv_filename = rdf_settings.get("path_or_buff", "rdf.csv")
        csv_file_path = os.path.join(working_dir, csv_filename)
        rdf_settings.update({"path_or_buff": csv_file_path})

        r_cut = rdf_settings.get("r_cut")

        bin_size = rdf_settings.get("bin_size")
        if isinstance(bin_size, (float, int)):
            bin_size = [bin_size]
        elif not isinstance(bin_size, (list, tuple)):
            pass
        mass = fw_spec.get("default_masses", [])
        if not mass:
            raise ValueError("Atomic masses not found")
        num_types = len(mass)
        num_mols = fw_spec.get("num_mols_list", [])
        num_atoms_per_mol = fw_spec.get("num_atoms_per_mol", [])
        partial_relations = rdf_settings.get("partial_relations", None)
        filename = rdf_settings.get("filename")
        save_mode = rdf_settings.get("save_mode")

        if rdf_type == "atomic":
            if use_default_atom_ids:
                num_mols = None
                num_atoms_per_mol = None
                if not partial_relations:
                    partial_relations = [[], []]
                    for i in range(1, num_types + 1):
                        for j in range(1, i + 1):
                            partial_relations[0].append(i)
                            partial_relations[1].append(j)
            else:
                mass = fw_spec.get("recalc_masses", [])
                if not num_mols or not num_atoms_per_mol:
                    raise ValueError(
                        "Number of molecules of each type and number of atoms "
                        "per molecule are not found")
                num_mols = [int(i) for i in num_mols]
                if not partial_relations:
                    partial_relations = [[], []]
                    for i in range(1, np.sum(num_atoms_per_mol) + 1):
                        for j in range(1, i + 1):
                            partial_relations[0].append(i)
                            partial_relations[1].append(j)
            atomic_working_dir = os.path.join(working_dir, rdf_type)
            os.makedirs(atomic_working_dir, exist_ok=True)
            for i, size in enumerate(bin_size):
                cur_working_dir = os.path.join(atomic_working_dir, str(size))
                os.makedirs(cur_working_dir, exist_ok=True)
                csv_file_path = os.path.join(cur_working_dir, csv_filename)
                calc_atomic_rdf(
                    r_cut,
                    size,
                    num_types,
                    mass,
                    partial_relations,
                    filename,
                    num_mols=num_mols,
                    num_atoms_per_mol=num_atoms_per_mol,
                    path_or_buff=csv_file_path,
                    save_mode=save_mode,
                )

        elif rdf_type == "molecular":
            if not num_mols or not num_atoms_per_mol:
                raise ValueError("Number of molecules of each type and number of atoms "
                                 "per molecule are not found")
            num_mols = [int(i) for i in num_mols]
            if not partial_relations:
                partial_relations = [[], []]
                for i in range(1, num_types + 1):
                    for j in range(1, len(num_atoms_per_mol) + 1):
                        partial_relations[0].append(i)
                        partial_relations[1].append(j)
            molecular_working_dir = os.path.join(working_dir, rdf_type)
            os.makedirs(molecular_working_dir, exist_ok=True)
            for i, size in enumerate(bin_size):
                cur_working_dir = os.path.join(molecular_working_dir, str(size))
                os.makedirs(cur_working_dir, exist_ok=True)
                csv_file_path = os.path.join(cur_working_dir, csv_filename)
                calc_molecular_rdf(
                    r_cut,
                    size,
                    num_types,
                    mass,
                    partial_relations,
                    filename,
                    num_mols=num_mols,
                    num_atoms_per_mol=num_atoms_per_mol,
                    path_or_buff=csv_file_path,
                    save_mode=save_mode,
                )

        rdf_settings_spec = {
            "r_cut": r_cut,
            "bin_size": bin_size,
            "num_types": num_types,
            "mass": mass,
            "partial_relations": partial_relations,
            "filename": filename,
            "num_mols": num_mols,
            "num_atoms_per_mol": num_atoms_per_mol,
            "path_or_buff": csv_filename,
            "save_mode": save_mode,
        }
        rdf_type_spec = rdf_type
        rdf_use_default_atom_ids_spec = use_default_atom_ids

        smiles_list = fw_spec.get("smiles", [])
        n_mols_dict = fw_spec.get("nmols", {})
        num_mols_list = num_mols
        num_atoms_per_mol_list = num_atoms_per_mol
        masses_list = mass
        lmp_box = fw_spec.get("box", None)

        return FWAction(
            update_spec={
                "cn_path": csv_file_path,
                "cn_settings": cn_settings_spec,
                "cn_type": cn_type_spec,
                "cn_use_default_atom_ids": cn_use_default_atom_ids_spec,
                "smiles": smiles_list,
                "nmols": n_mols_dict,
                "num_mols_list": num_mols_list,
                "num_atoms_per_mol": num_atoms_per_mol_list,
                "masses": masses_list,
                "box": lmp_box,
            }
        )
