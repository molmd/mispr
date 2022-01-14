# coding: utf-8


# Defines firetasks for reading output of Ambertools programs and LAMMPS.

import os
import json
import inspect

import numpy as np
import pandas as pd
import logging
from scipy.signal import argrelmin, find_peaks

from fireworks.utilities.fw_serializers import DATETIME_HANDLER
from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from pymatgen.io.ambertools import PrmtopParser
from pymatgen.core.structure import Molecule

from mdproptools.structural.rdf_cn import (
    calc_atomic_rdf,
    calc_molecular_rdf,
    calc_atomic_cn,
    calc_molecular_cn,
)
from mdproptools.structural.cluster_analysis import (
    get_clusters,
    get_unique_configurations,
)
from mdproptools.dynamical.diffusion import Diffusion

from mispr.lammps.utilities.utilities import (
    get_db,
    process_ff_doc,
    add_ff_labels_to_dict,
)

from mispr.lammps.defaults import (
    MSD_SETTINGS,
    DIFF_SETTINGS,
    RDF_SETTINGS,
    CN_SETTINGS,
    CLUSTERS_SETTINGS,
)
from mispr.gaussian.utilities.metadata import get_mol_formula

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

GAFF_DOI = "https://doi.org/10.1002/jcc.20035"


@explicit_serialize
class ProcessPrmtop(FiretaskBase):
    _fw_name = "Process Prmtop"
    required_params = []
    optional_params = [
        "molecule" "working_dir",
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
        save_to_db = self.get("save_ff_to_db", False)
        save_to_file = self.get("save_ff_to_file", True)

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

        molecule = self.get("molecule", fw_spec["prev_calc_molecule"])

        if not isinstance(molecule, Molecule):
            raise Exception('"molecule" is not a pymatgen Molecule object')

        # TODO: decide what unique id for the ff labels. If staying as smiles \
        #  string, decide how to get smiles.
        unique_mol_name = self.get("unique_molecule_name", get_mol_formula(molecule))
        ff_param_dict_general = PrmtopParser(prmtop_file_path, molecule, "").to_dict()
        ff_param_dict_system = add_ff_labels_to_dict(
            ff_param_dict_general, unique_mol_name
        )

        gaff_doi = self.get("gaff_doi", GAFF_DOI)

        if save_to_db:
            ff_doc = ff_param_dict_general.copy()
            ff_db = get_db(input_db=db)
            ff_db.insert_force_field(ff_doc, "gaff", doi=gaff_doi)

        if save_to_file:
            ff_doc = process_ff_doc(ff_param_dict_general, "gaff", doi=gaff_doi)
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
            propagate=True,
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
        outputs_dir = diff_settings.get(
            "outputs_dir", os.path.abspath(os.path.join(working_dir, "..", "..", "nvt"))
        )
        msd_df = None
        diff = Diffusion(
            timestep=diff_settings["timestep"],
            units=diff_settings["units"],
            outputs_dir=outputs_dir,
            diff_dir=working_dir,
        )
        if msd_method == "from_dump":
            num_mols = diff_settings.pop("num_mols", fw_spec.get("num_mols_list", []))
            num_atoms_per_mol = diff_settings.pop("num_atoms_per_mol", fw_spec.get("num_atoms_per_mol", []))
            num_mols = [int(i) for i in num_mols] if num_mols else None
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
                },
            )
        elif msd_method == "from_log":
            file_pattern = diff_settings.get("file_pattern", "log.lammps*")
            msd_df = diff.get_msd_from_log(file_pattern)

        else:
            raise ValueError(
                f"Unsupported msd calculation method: {msd_method}"
                f"Supported methods are from_dump and from_log"
            )

        diffusion = diff.calc_diff(
            msd_df[0],
            **{
                i: j
                for i, j in diff_settings.items()
                if i in inspect.getfullargspec(diff.calc_diff).args
            },
        )

        diff_list = diffusion["diffusion (m2/s)"].to_list()
        std_list = diffusion["std"].to_list()
        r2_list = diffusion["R2"].to_list()

        if (
            msd_method == "from_dump"
            and diff_settings["avg_interval"]
            and diff_settings["diff_dist"]
        ):
            diff.get_diff_dist(
                msd_df[2],
                **{
                    i: j
                    for i, j in diff_settings.items()
                    if i in inspect.getfullargspec(diff.get_diff_dist).args
                },
            )

        diff_settings.update(
            {
                "diff_path": working_dir,
                "diffusion": diff_list,
                "std": std_list,
                "r2": r2_list,
            }
        )
        return FWAction(
            update_spec={
                "smiles": fw_spec.get("smiles", []),
                "nmols": fw_spec.get("nmols", {}),
                "num_mols_list": fw_spec.get("num_mols_list", []),
                "box": fw_spec.get("box", None),
                "num_atoms_per_mol": fw_spec.get("num_atoms_per_mol", []),
                "default_masses": fw_spec.get("default_masses", []),
                "recalc_masses": fw_spec.get("recalc_masses", []),
                "diffusion": diff_settings,
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
        mass = rdf_settings.get("mass", fw_spec.get("default_masses", []))
        if not mass:
            raise ValueError("Atomic masses not found")
        num_types = len(mass)
        num_mols = rdf_settings.get("num_mols", fw_spec.get("num_mols_list", []))
        num_atoms_per_mol = rdf_settings.get("num_atoms_per_mol", fw_spec.get("num_atoms_per_mol", []))
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
                        "per molecule are not found"
                    )
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
                raise ValueError(
                    "Number of molecules of each type and number of atoms "
                    "per molecule are not found"
                )
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
            "rdf_type": rdf_type,
            "rdf_use_default_atom_ids": use_default_atom_ids,
            "rdf_path": csv_file_path,
        }

        return FWAction(
            update_spec={
                "smiles": fw_spec.get("smiles", []),
                "nmols": fw_spec.get("nmols", {}),
                "num_mols_list": num_mols,
                "num_atoms_per_mol": num_atoms_per_mol,
                "masses": mass,
                "box": fw_spec.get("box", None),
                "rdf": rdf_settings_spec,
            }
        )


@explicit_serialize
class CalcCN(FiretaskBase):
    _fw_name = "Calculate CN"
    required_params = []
    optional_params = ["cn_settings", "working_dir"]

    def run_task(self, fw_spec):
        # Get inputs to cn function from defaults or from user, which then
        # updates the defaults
        cn = None
        cn_settings = CN_SETTINGS.copy()
        cn_settings.update(self.get("cn_settings", {}))

        # atomic or molecular (for CoM)
        cn_type = cn_settings.get("cn_type", "atomic")

        working_dir = self.get("working_dir", os.getcwd())
        os.makedirs(working_dir, exist_ok=True)

        use_default_atom_ids = cn_settings.get("use_default_atom_ids", False)

        csv_filename = cn_settings.get("path_or_buff", "cn.csv")
        csv_file_path = os.path.join(working_dir, csv_filename)
        cn_settings.update({"path_or_buff": csv_file_path})

        # check if user provided the cutoff radius to use in calculating the CN
        r_cut = cn_settings.get("r_cut", None)
        # if the cutoff radius is not provided, attempt to automatically detect it from
        # the rdf
        if not r_cut:
            rdf_path = fw_spec.get("rdf", {}).get("rdf_path")
            rdf_type = fw_spec.get("rdf", {}).get("rdf_type")
            rdf_use_default_atom_ids = fw_spec.get("rdf", {}).get(
                "rdf_use_default_atom_ids"
            )
            if not rdf_path:
                raise ValueError(
                    "Cutoff distance required for calculating the CN is "
                    "not found and cannot be computed due to a missing "
                    "RDF file"
                )
            else:
                if (
                    cn_type != rdf_type
                    or use_default_atom_ids != rdf_use_default_atom_ids
                ):
                    raise ValueError("CN settings do not match those of RDF")
                else:
                    r_cut = []
                    rdf_df = pd.read_csv(rdf_path)
                    for col in rdf_df:
                        if col != "r ($\AA$)" and col != "g_full(r)":
                            new_df = rdf_df[["r ($\AA$)", col]]
                            new_df.columns = ["x", "y"]
                            peaks = find_peaks(new_df["y"], height=1)[0]
                            minima = argrelmin(new_df["y"].values)[0]
                            min_3 = new_df.loc[minima[minima > peaks[0]][0:3]]
                            r_cut.append(np.mean(min_3["x"].tolist()))

        bin_size = cn_settings.get("bin_size")
        if isinstance(bin_size, (float, int)):
            bin_size = [bin_size]
        elif not isinstance(bin_size, (list, tuple)):
            pass
        mass = cn_settings.get("mass", fw_spec.get("default_masses", []))
        if not mass:
            raise ValueError("Atomic masses not found")
        num_types = len(mass)
        num_mols = cn_settings.get("num_mols", fw_spec.get("num_mols_list", []))
        num_atoms_per_mol = cn_settings.get("num_atoms_per_mol", fw_spec.get("num_atoms_per_mol", []))
        partial_relations = cn_settings.get("partial_relations", None)
        filename = cn_settings.get("filename")
        save_mode = cn_settings.get("save_mode")

        if cn_type == "atomic":
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
                        "per molecule are not found"
                    )
                num_mols = [int(i) for i in num_mols]
                if not partial_relations:
                    partial_relations = [[], []]
                    for i in range(1, np.sum(num_atoms_per_mol) + 1):
                        for j in range(1, i + 1):
                            partial_relations[0].append(i)
                            partial_relations[1].append(j)
            atomic_working_dir = os.path.join(working_dir, cn_type)
            os.makedirs(atomic_working_dir, exist_ok=True)
            for i, size in enumerate(bin_size):
                cur_working_dir = os.path.join(atomic_working_dir, str(size))
                os.makedirs(cur_working_dir, exist_ok=True)
                csv_file_path = os.path.join(cur_working_dir, csv_filename)
                cn = calc_atomic_cn(
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

        elif cn_type == "molecular":
            if not num_mols or not num_atoms_per_mol:
                raise ValueError(
                    "Number of molecules of each type and number of atoms "
                    "per molecule are not found"
                )
            num_mols = [int(i) for i in num_mols]
            if not partial_relations:
                partial_relations = [[], []]
                for i in range(1, num_types + 1):
                    for j in range(1, len(num_atoms_per_mol) + 1):
                        partial_relations[0].append(i)
                        partial_relations[1].append(j)
            molecular_working_dir = os.path.join(working_dir, cn_type)
            os.makedirs(molecular_working_dir, exist_ok=True)
            for i, size in enumerate(bin_size):
                cur_working_dir = os.path.join(molecular_working_dir, str(size))
                os.makedirs(cur_working_dir, exist_ok=True)
                csv_file_path = os.path.join(cur_working_dir, csv_filename)
                cn = calc_molecular_cn(
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

        cn_settings_spec = {
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
            "cn_type": cn_type,
            "cn_use_default_atom_ids": use_default_atom_ids,
            "cn_path": csv_file_path,
            "cn": cn.loc[0].to_dict(),
        }

        return FWAction(
            update_spec={
                "smiles": fw_spec.get("smiles", []),
                "nmols": fw_spec.get("nmols", {}),
                "num_mols_list": num_mols,
                "num_atoms_per_mol": num_atoms_per_mol,
                "masses": mass,
                "box": fw_spec.get("box", None),
                "cn": cn_settings_spec,
            }
        )


@explicit_serialize
class ExtractClusters(FiretaskBase):
    _fw_name = "Extract Atomic Clusters"
    required_params = []
    optional_params = ["cluster_settings", "working_dir"]

    def run_task(self, fw_spec):
        cluster_settings = CLUSTERS_SETTINGS.copy()
        cluster_settings.update(self.get("cluster_settings", {}))
        working_dir = self.get("working_dir", os.getcwd())
        os.makedirs(working_dir, exist_ok=True)
        filename = cluster_settings.get("filename")
        atom_type = cluster_settings.get("atom_type")
        if not atom_type:
            raise ValueError("No atom type specified to perform cluster analysis")
        num_mols = cluster_settings.get("num_mols", fw_spec.get("num_mols_list", []))
        num_mols = [int(i) for i in num_mols]
        num_atoms_per_mol = cluster_settings.get(
            "num_atoms_per_mol", fw_spec.get("num_atoms_per_mol", [])
        )
        # TODO: get elements from fw_spec
        elements = cluster_settings.get("elements", None)
        prev_rcut = fw_spec.get("cn", {}).get("r_cut")
        cur_rcut = min(prev_rcut) if prev_rcut else None
        r_cut = cluster_settings.get("r_cut", cur_rcut)
        full_trajectory = cluster_settings.get("full_trajectory", True)
        frame = cluster_settings.get("frame", None)
        alter_atom_ids = cluster_settings.get("alter_atom_ids", False)
        max_force = cluster_settings.get("max_force", 0.75)
        cluster_count = get_clusters(
            filename=filename,
            atom_type=atom_type,
            r_cut=r_cut,
            full_trajectory=full_trajectory,
            frame=frame,
            num_mols=num_mols,
            num_atoms_per_mol=num_atoms_per_mol,
            elements=elements,
            alter_atom_ids=alter_atom_ids,
            max_force=max_force,
            working_dir=working_dir,
        )

        # TODO: if molecules are taken from fw_spec, need to make sure if mol_names
        #  are provided they match the order in fw_spec (which is sorted by the
        #  num_atoms/mol)
        type_coord_atoms = cluster_settings.get("type_coord_atoms", None)
        perc = cluster_settings.get("perc", None)
        cum_perc = cluster_settings.get("cum_perc", 90)
        clusters, configurations = get_unique_configurations(
            cluster_pattern="Cluster_*",
            r_cut=r_cut,
            molecules=cluster_settings.get("molecules", fw_spec.get("molecules")),
            type_coord_atoms=type_coord_atoms,
            working_dir=working_dir,
            find_top=cluster_settings.get("find_top", True),
            perc=perc,
            cum_perc=cum_perc,
            mol_names=cluster_settings.get("mol_names", None),
            zip=cluster_settings.get("zip", True),
        )
        top_config_files = [
            f"{working_dir}/{file}"
            for file in os.listdir(working_dir)
            if file.startswith("conf_") and file.endswith(".xyz")
        ]

        cluster_analysis_spec = {
            "atom_type": atom_type,
            "r_cut": r_cut,
            "full_trajectory": full_trajectory,
            "frame": frame,
            "num_mols": num_mols,
            "num_atoms_per_mol": num_atoms_per_mol,
            "elements": elements,
            "alter_atom_ids": alter_atom_ids,
            "max_force": max_force,
            "num_clusters": cluster_count,
            "num_configurations": len(configurations),
            "type_coord_atoms": type_coord_atoms,
            "perc": perc,
            "cum_perc": cum_perc,
            "filename": filename,
            "working_dir": working_dir,
        }

        return FWAction(
            update_spec={
                "clusters": cluster_analysis_spec,
                "top_config_files": top_config_files,
            }
        )


@explicit_serialize
class ProcessAnalysis(FiretaskBase):
    _fw_name = "Process Analysis Calculations"
    required_params = ["analysis_list"]
    optional_params = [
        "db",
        "save_analysis_to_db",
        "save_analysis_to_file",
        "working_dir",
    ]

    def run_task(self, fw_spec):
        db = self.get("db", None)
        working_dir = self.get("working_dir", fw_spec.get("working_dir", os.getcwd()))
        save_analysis_to_db = self.get("save_analysis_to_db", False)
        save_analysis_to_file = self.get("save_analysis_to_file", True)

        box = fw_spec.get("box", None)
        if box:
            box = box.as_dict()

        systems_dict = {
            "smiles": fw_spec.get("smiles", []),
            "nmols": fw_spec.get("num_mols_list", []),
            "natoms_per_mol": fw_spec.get("num_atoms_per_mol", []),
            "mass": fw_spec.get("masses", []),
            "box": box,
        }

        for analysis in self["analysis_list"]:
            systems_dict[analysis] = fw_spec.get(analysis)

        if fw_spec.get("run_id_list"):
            systems_dict["run_ids"] = fw_spec["run_id_list"]

        if save_analysis_to_db:
            db = get_db(input_db=db)
            db.insert_system(systems_dict)

        if fw_spec.get("lammps_run_loc_list"):
            systems_dict["lammps_run_locs"] = fw_spec["lammps_run_loc_list"]

        if save_analysis_to_file:
            if "lammps_run_ids" in systems_dict:
                del systems_dict["lammps_run_ids"]
            file = os.path.join(working_dir, "system.json")
            with open(file, "w") as f:
                f.write(json.dumps(systems_dict, default=DATETIME_HANDLER))

        logger.info("Analysis calculations complete")
