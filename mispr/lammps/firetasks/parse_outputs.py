# coding: utf-8

# Defines firetasks for reading output of Ambertools programs and LAMMPS

import os
import json
import inspect
import numpy as np
import pandas as pd

import pymatgen.core.structure as pmgcs
import pymatgen.io.ambertools as pmgia
import infrastructure.gaussian.utils.utils as iguu
import infrastructure.lammps.utils.utils as iluu
import analysis.lammps.dynamical.diffusion as aldd
import analysis.lammps.structural.rdf_cn as alsrc
import fireworks.core.firework as fwcfw
import fireworks.utilities.fw_utilities as fwuu
from fireworks.user_objects.firetasks.script_task import PyTask

__author__ = 'Matthew Bliss'
__version__ = '0.0.1'
__email__ = 'matthew.bliss@tufts.edu'
__date__ = 'Apr 15, 2020'

GAFF_DOI = "https://doi.org/10.1002/jcc.20035"

DEFAULT_RDF_SETTINGS = {"r_cut": 20,
                        "bin_size": [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5],
                        "filename": os.path.abspath(
                            os.path.join("../../../lammps", "nvt", "log.lammps")
                        ),
                        "path_or_buff": "rdf.csv",
                        "save_mode": True}

DEFAULT_MSD_SETTINGS = {"dt": 1,
                        "tao_coeff": 4,
                        "msd_type": "com",
                        "return_all": False,
                        "com_drift": False,
                        "avg_interval": False,
                        "save_msd": True}

DEFAULT_DIFF_SETTINGS = {"initial_time": None,
                         "final_time": None,
                         "dimension": 3}


@fwuu.explicit_serialize
class ProcessPrmtop(fwcfw.FiretaskBase):
    _fw_name = "Process Prmtop"
    required_params = ["molecule"]
    optional_params = ["working_dir", "db", "prmtop_path", "prmtop_filename",
                       "prmtop_dir", "unique_molecule_name",
                       "system_force_field_dict", "doi", "save_ff_to_db",
                       "save_ff_to_file", "ff_filename"]

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
            raise Exception('Neither name nor path to a prmtop file was \
                            provided; either provide prmtop file name as \
                            "prmtop_filename" or path to prmtop as \
                            "prmtop_path".')

        if not os.path.exists(prmtop_file_path):
            raise Exception('"prmtop_file_path" is not a valid path; check \
                            that "prmtop_path" is correct or that \
                            "prmtop_filename" exists in "prmtop_dir"')

        molecule = self.get("molecule")

        if not isinstance(molecule, pmgcs.Molecule):
            raise Exception('"molecule" is not a pymatgen Molecule object')

        # TODO: decide what unique id for the ff labels. If staying as smiles \
        #  string, decide how to get smiles.
        unique_mol_name = self.get("unique_molecule_name",
                                   iguu.get_mol_formula(molecule))
        # print("UNIQUE MOL NAME")
        # print(unique_mol_name)
        ff_param_dict_general = pmgia.PrmtopParser(prmtop_file_path,
                                                   molecule, '').to_dict()
        ff_param_dict_system = iluu.add_ff_labels_to_dict(
            ff_param_dict_general, unique_mol_name
        )

        gaff_doi = self.get("gaff_doi", GAFF_DOI)

        if db and save_to_db:
            ff_db = iluu.get_db(input_db=db)
            ff_db.insert_force_field(ff_param_dict_general, "gaff",
                                     doi=gaff_doi)

        if save_to_file:
            ff_doc = iluu.process_ff_doc(ff_param_dict_general, "gaff",
                                         doi=gaff_doi)
            with open(os.path.join(working_dir, ff_filename), "w") as file:
                json.dump(ff_doc, file)

        sys_ff_dict = fw_spec.get("system_force_field_dict",
                                  self.get("system_force_field_dict", {}))
        sys_ff_dict[unique_mol_name] = ff_param_dict_system

        return fwcfw.FWAction(
            stored_data={'ff_param_dict_system': ff_param_dict_system,
                         'ff_param_dict_general': ff_param_dict_general},
            update_spec={"system_force_field_dict": sys_ff_dict}
        )
                        # mod_spec = [{'_set': {"system_force_field_dict": \
                        #                                       sys_ff_dict}}])


@fwuu.explicit_serialize
class GetMSD(fwcfw.FiretaskBase):
    _fw_name = "Get MSD"
    required_params = ["msd_method"]
    optional_params = ["working_dir", "msd_settings"]

    def run_task(self, fw_spec):
        method = self.get("msd_method")
        working_dir = self.get("working_dir", os.getcwd())
        os.makedirs(working_dir, exist_ok=True)

        msd_settings = DEFAULT_MSD_SETTINGS.copy()
        msd_settings.update(self.get("msd_settings", {}))

        if method == "from_dump":
            mass = fw_spec.get("default_masses", None)
            num_mols = fw_spec.get("num_mols_list", None)
            num_atoms_per_mol = fw_spec.get("num_atoms_per_mol", None)
            mass = fw_spec.get("default_masses", None)

            file_pattern = msd_settings.get("file_pattern",
                                            os.path.abspath(os.path.join(working_dir,
                                                                         "../../../lammps", "nvt", "dump.nvt.*.dump")))

            msd = aldd.get_msd_from_dump(
                file_pattern,
                num_mols=num_mols,
                num_atoms_per_mol=num_atoms_per_mol,
                mass=mass,
                working_dir=working_dir,
                **{i: j for i, j in msd_settings.items() if i in
                   inspect.getfullargspec(aldd.get_msd_from_dump).args}
            )
        elif method == "from_log":
            file_pattern = msd_settings.get("file_pattern",
                                    os.path.abspath(
                                        os.path.join(working_dir, "../../../lammps", "nvt",
                                                     "log.lammps")))
            print(file_pattern)
            msd = aldd.get_msd_from_log(
                file_pattern,
                working_dir=working_dir,
                **{i: j for i, j in msd_settings.items() if i in
                   inspect.getfullargspec(aldd.get_msd_from_log).args}
            )

        smiles_list = fw_spec.get("smiles", [])
        n_mols_dict = fw_spec.get("nmols", {})
        num_mols_list = fw_spec.get("num_mols_list", [])
        lmp_box = fw_spec.get("box", None)
        num_atoms_per_mol_list = fw_spec.get("num_atoms_per_mol", [])
        default_masses_list = fw_spec.get("default_masses", [])
        recalc_masses_list = fw_spec.get("recalc_masses", [])

        return fwcfw.FWAction(
            update_spec={"msd_file_path": os.path.join(working_dir, "msd.csv"),
                         "smiles": smiles_list,
                         "nmols": n_mols_dict,
                         "num_mols_list": num_mols_list,
                         "box": lmp_box,
                         "num_atoms_per_mol":
                             num_atoms_per_mol_list,
                         "default_masses": default_masses_list,
                         "recalc_masses": recalc_masses_list}
        )


@fwuu.explicit_serialize
class CalcDiff(fwcfw.FiretaskBase):
    _fw_name = "Calculate Diffusion"
    required_params = []
    optional_params = ["msd", "working_dir", "diff_settings", "msd_file_path"]

    def run_task(self, fw_spec):
        working_dir = self.get("working_dir", os.getcwd())
        os.makedirs(working_dir, exist_ok=True)

        msd = self.get("msd")
        msd_file_path = fw_spec.get(
            "msd_file_path",
            self.get("msd_file_path",
                     os.path.abspath(os.path.join(working_dir,
                                                  "../../../lammps", "msd"
                                                  "msd.csv")))
        )

        if not msd:
            msd = pd.read_csv(msd_file_path)

        diff_settings = DEFAULT_DIFF_SETTINGS.copy()
        diff_settings.update(self.get("diff_settings", {}))
        diff_settings.update({"working_dir": working_dir})

        diff, models = aldd.get_diff(
            msd, **{i: j for i, j in diff_settings.items() if i in
                    inspect.getfullargspec(aldd.get_diff).args}
        )

        smiles_list = fw_spec.get("smiles", [])
        n_mols_dict = fw_spec.get("nmols", {})
        num_mols_list = fw_spec.get("num_mols_list", [])
        lmp_box = fw_spec.get("box", None)
        num_atoms_per_mol_list = fw_spec.get("num_atoms_per_mol", [])
        default_masses_list = fw_spec.get("default_masses", [])
        recalc_masses_list = fw_spec.get("recalc_masses", [])

        return fwcfw.FWAction(update_spec={"diffusion_calc_dir": working_dir,
                                           "smiles": smiles_list,
                                           "nmols": n_mols_dict,
                                           "num_mols_list": num_mols_list,
                                           "box": lmp_box,
                                           "num_atoms_per_mol":
                                               num_atoms_per_mol_list,
                                           "default_masses":
                                               default_masses_list,
                                           "recalc_masses":
                                               recalc_masses_list})


@fwuu.explicit_serialize
class GetRDF(fwcfw.FiretaskBase):
    _fw_name = "Get RDF"
    required_params = []
    optional_params = ["rdf_type", "rdf_settings", "working_dir",
                       "use_default_atom_ids"]

    def run_task(self, fw_spec):

        # Get inputs to rdf function from defaults or from user, which then
        # updates the defaults
        init_rdf_settings = self.get("rdf_settings", DEFAULT_RDF_SETTINGS)
        print(init_rdf_settings)
        rdf_settings = DEFAULT_RDF_SETTINGS.copy()
        rdf_settings.update(init_rdf_settings)

        # atomic or molecular (for CoM)
        rdf_type = rdf_settings.get("rdf_type", self.get("rdf_type", "atomic"))

        working_dir = self.get("working_dir", os.getcwd())
        os.makedirs(working_dir, exist_ok=True)

        print(working_dir)

        use_default_atom_ids = rdf_settings.get(
            "use_default_atom_ids",
            self.get("use_default_atom_ids",
                     False
        ))

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
        print(mass)
        if not mass:
            # TODO: add exception
            pass

        num_types = len(mass)

        partial_relations = rdf_settings.get("partial_relations", None)

        filename = rdf_settings.get("filename")

        if use_default_atom_ids:
            mass = fw_spec.get("default_masses", [])
            num_mols = None
            num_atoms_per_mol = None

            if not mass:
                pass

            if partial_relations is None and rdf_type == "atomic":
                partial_relations = [[], []]
                for i in range(1, num_types + 1):
                    if rdf_type == "atomic":
                        for j in range(1, i + 1):
                            partial_relations[0].append(i)
                            partial_relations[1].append(j)
                    elif rdf_type == "molecular":
                        num_species = len(fw_spec.get("smiles"))
                        for j in range(1, num_species + 1):
                            partial_relations[0].append(i)
                            partial_relations[1].append(j)
        else:
            mass = fw_spec.get("recalc_masses", [])
            num_mols = fw_spec.get("num_mols_list", [])
            num_atoms_per_mol = fw_spec.get("num_atoms_per_mol", [])

            if not num_mols or not num_atoms_per_mol:
                # TODO: add exception
                pass

            if partial_relations is None:
                partial_relations = [[], []]
                for i in range(1, np.sum(num_atoms_per_mol) + 1):
                    if rdf_type == "atomic":
                        for j in range(1, i + 1):
                            partial_relations[0].append(i)
                            partial_relations[1].append(j)
                    elif rdf_type == "molecular":
                        for j in range(1, len(num_atoms_per_mol) + 1):
                            partial_relations[0].append(i)
                            partial_relations[1].append(j)

        save_mode = rdf_settings.get("save_mode")

        print(partial_relations)

        if rdf_type == "atomic":
            atomic_working_dir = os.path.join(working_dir, rdf_type)
            os.makedirs(atomic_working_dir, exist_ok=True)
            for i, size in enumerate(bin_size):
                cur_working_dir = os.path.join(atomic_working_dir, str(size))
                os.makedirs(cur_working_dir, exist_ok=True)
                csv_file_path = os.path.join(cur_working_dir, csv_filename)
                rdf_data = alsrc.calc_atomic_rdf(r_cut,
                                                 size,
                                                 num_types,
                                                 mass,
                                                 partial_relations,
                                                 filename,
                                                 num_mols=num_mols,
                                                 num_atoms_per_mol=
                                                            num_atoms_per_mol,
                                                 path_or_buff=csv_file_path,
                                                 save_mode=save_mode)
        elif rdf_type == "molecular":
            molecular_working_dir = os.path.join(working_dir, rdf_type)
            os.makedirs(molecular_working_dir, exist_ok=True)
            for i, size in enumerate(bin_size):
                cur_working_dir = os.path.join(molecular_working_dir,
                                               str(size))
                os.makedirs(cur_working_dir, exist_ok=True)
                csv_file_path = os.path.join(cur_working_dir, csv_filename)
                rdf_data = alsrc.calc_molecular_rdf(r_cut,
                                                    size,
                                                    num_types,
                                                    mass,
                                                    partial_relations,
                                                    filename,
                                                    num_mols=num_mols,
                                                    num_atoms_per_mol=
                                                            num_atoms_per_mol,
                                                    path_or_buff=csv_file_path,
                                                    save_mode=save_mode)
        rdf_settings_spec = {"r_cut": r_cut,
                             "bin_size": bin_size,
                             "num_types": num_types,
                             "mass": mass,
                             "partial_relations": partial_relations,
                             "filename": filename,
                             "num_mols": num_mols,
                             "num_atoms_per_mol": num_atoms_per_mol,
                             "path_or_buff": csv_filename,
                             "save_mode": save_mode}
        rdf_type_spec = rdf_type
        rdf_use_default_atom_ids_spec = use_default_atom_ids

        smiles_list = fw_spec.get("smiles", [])
        n_mols_dict = fw_spec.get("nmols", {})
        num_mols_list = num_mols
        num_atoms_per_mol_list = num_atoms_per_mol
        masses_list = mass
        lmp_box = fw_spec.get("box", None)

        return fwcfw.FWAction(update_spec={
            "rdf_calc_dir": working_dir,
            "rdf_settings": rdf_settings_spec,
            "rdf_type": rdf_type_spec,
            "rdf_use_default_atom_ids": rdf_use_default_atom_ids_spec,
            "smiles": smiles_list,
            "nmols": n_mols_dict,
            "num_mols_list": num_mols_list,
            "num_atoms_per_mol":
                num_atoms_per_mol_list,
            "masses": masses_list,
            "box": lmp_box
        })


if __name__ == "__main__":
    # test_string = "test_string\t\n"
    # print([test_string])
    # print([test_string.strip()])
    # args_dict = {"first": 1, "second": 2}
    # for k, v in args_dict.items():
    #     print(k, v)
    # diff_dir_path = "/Users/matt/Documents/Research/DHPS/Depro_DHPS/Lammps/SPCE/Viscosity/0.1_M"
    # test_msd = aldd.get_msd_from_log("log.lammps", dt = 2, working_dir = diff_dir_path)
    # print(test_msd.head(10))
    # print(test_msd.tail(10))
    # diff, model = aldd.get_diff(test_msd, working_dir = diff_dir_path)
    # print(diff)

    # params = inspect.getfullargspec(aldd.get_msd_from_dump)
    # print(type(params))
    # print(params)
    # print(params.args)
    # print(type(params.args))
    # print(params.defaults)

    string = "log.lammps"
    forslash = "\\"
    print(string)
    print(forslash)
    if forslash not in string and "/" not in string:
        print(False)

    # os.chdir("/Users/matt/Documents/Research/Infrastructure/case_study")
    # with open("test_2.txt", "w") as f:
    #     f.write("hello")
    test_list = [1.0, 2.1]
    print(type(test_list))
    print(isinstance(test_list, list))