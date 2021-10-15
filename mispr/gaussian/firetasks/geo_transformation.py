# coding: utf-8


# Defines firetasks for performing various molecule transformations.

import os
import copy
import logging
import itertools

from fireworks import Workflow
from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.analysis.fragmenter import open_ring

from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.utilities.metadata import get_mol_formula
from mispr.gaussian.utilities.db_utilities import get_db

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

DEFAULT_KEY = "gout_key"


@explicit_serialize
class ProcessMoleculeInput(FiretaskBase):
    required_params = ["mol"]
    optional_params = [
        "operation_type",
        "db",
        "save_to_db",
        "charge",
        "update_duplicates",
        "save_to_file",
        "fmt",
        "filename",
        "from_fw_spec",
        "local_opt",
        "force_field",
        "steps",
        "str_type",
        "working_dir"
    ]

    @staticmethod
    def _run_to_mol_object(run):
        return Molecule.from_dict(run["output"]["output"]["molecule"])

    @staticmethod
    def _from_fw_spec(mol, fw_spec):
        # mol = key in this case
        available_runs = fw_spec["gaussian_output"]
        if not isinstance(mol, dict):
            mol = available_runs.get(mol, mol)
        else:
            if isinstance(mol["mol"], list):
                mol["mol"] = [
                    ProcessMoleculeInput._from_fw_spec(i, fw_spec) for i in mol["mol"]
                ]
            else:
                mol["mol"] = ProcessMoleculeInput._from_fw_spec(mol["mol"], fw_spec)
        return mol

    def run_task(self, fw_spec):
        mol = self["mol"]
        operation_type = self.get("operation_type", "get_from_mol")
        working_dir = self.get("working_dir") or os.getcwd()
        db = self.get("db")

        if self.get("from_fw_spec"):
            mol = self._from_fw_spec(mol, fw_spec)

        output_mol = process_mol(
            operation_type=operation_type,
            mol=mol,
            working_dir=working_dir,
            db=db,
            local_opt=self.get("local_opt", False),
            force_field=self.get("force_field"),
            steps=self.get("steps"),
            charge=self.get("charge"),
            str_type=self.get("str_type"),
        )

        if self.get("save_to_db"):
            db = get_db(db) if db else get_db()
            update_duplicates = self.get("update_duplicates", False)
            db.insert_molecule(output_mol, update_duplicates=update_duplicates)

        if self.get("save_to_file"):
            fmt = self.get("fmt", "xyz")
            filename = self.get("filename", get_mol_formula(output_mol))
            file = os.path.join(working_dir, f"{filename}.{fmt}")
            output_mol.to(fmt, file)
        fw_spec["prev_calc_molecule"] = output_mol  # Note: This should ideally
        # be part of FWaction, however because mpi doesn"t support pymatgen, we
        # should be careful about what is being passed to the next firework


@explicit_serialize
class ConvertToMoleculeObject(FiretaskBase):
    """
    Reads a molecule from a file, converts it to a mol object,
    and saves it as dict to mongodb.
    Supported file formats include
        xyz|pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv,
        gaussian input (gjf|g03|g09|com|inp),
        Gaussian output (.out), and
        pymatgen"s JSON serialized molecules.
    Requires openbabel to be installed.
    """

    required_params = ["mol_file"]
    optional_params = ["db", "save_to_db", "update_duplicates"]

    def run_task(self, fw_spec):
        working_dir = os.getcwd()
        file_name = self["mol_file"]
        mol = process_mol(file_name, working_dir)
        if self.get("save_to_db", True):
            mol_db = get_db(self.get("db"))
            mol_db.insert_molecule(
                mol, update_duplicates=self.get("update_duplicates", False)
            )
        fw_spec["prev_calc_molecule"] = mol  # Note: This should ideally be
        # part of FWaction, however because mpi doesn"t support pymatgen, we
        # should be careful about what is being passed to the next firework


@explicit_serialize
class RetrieveMoleculeObject(FiretaskBase):
    """
    Returns a molecule object from the database using inchi as an identifier
    """

    required_params = ["inchi"]
    optional_params = ["db", "save_mol_file", "fmt", "filename"]

    def run_task(self, fw_spec):
        inchi = self["inchi"]
        mol = process_mol("get_from_mol_db", inchi, db=self.get("db"))
        if mol and self.get("save_mol_file", False):
            working_dir = os.getcwd()
            file_name = self.get(
                "filename.{}".format(self.get("fmt", "xyz")),
                "{}.{}".format(get_mol_formula(mol), self.get("fmt", "xyz")),
            )
            mol_file = os.path.join(working_dir, file_name)
            mol.to(self.get("fmt", "xyz"), mol_file)
        fw_spec["prev_calc_molecule"] = mol


@explicit_serialize
class AttachFunctionalGroup(FiretaskBase):
    """
    Attaches a functional group to a molecule; requires the name of the functional
    group and a molecule (from a prev. calc or as a molecule object)
    """

    required_params = ["func_grp", "index"]
    optional_params = [
        "db",
        "molecule",
        "bond_order",
        "save_to_db",
        "update_duplicates",
        "save_mol_file",
        "fmt",
        "filename",
    ]

    def run_task(self, fw_spec):
        db = get_db(self.get("db"))
        if fw_spec.get("prev_calc_molecule"):
            mol = fw_spec.get("prev_calc_molecule")
        elif self.get("molecule"):
            mol = self.get("molecule")
        func_grp_dict = db.retrieve_fg(self["func_grp"])
        func_grp = Molecule(func_grp_dict["species"], func_grp_dict["coords"])
        derived_mol = Molecule.copy(mol)
        derived_mol.substitute(index=self["index"], func_grp=func_grp)
        if self.get("save_to_db", True):
            db.insert_derived_mol(
                derived_mol, update_duplicates=self.get("update_duplicates", False)
            )
        if self.get("save_mol_file", False):
            working_dir = os.getcwd()
            file_name = self.get("filename", get_mol_formula(derived_mol))
            file_name = ("{}.{}".format(file_name, self.get("fmt", "xyz")),)
            derived_mol_file = os.path.join(working_dir, file_name)
            derived_mol.to(self.get("fmt", "xyz"), derived_mol_file)
        fw_spec["prev_calc_molecule"] = derived_mol


@explicit_serialize
class LinkMolecules(FiretaskBase):
    """
    Links two molecules using one site from the first and another site from the
    second molecule.
    """

    required_params = ["index1", "index2"]
    optional_params = [
        "db",
        "mol1",
        "mol2",
        "bond_order",
        "save_to_db",
        "update_duplicates",
        "save_mol_file",
        "fmt",
        "filename",
    ]

    def run_task(self, fw_spec):
        db = get_db(self.get("db"))
        mol1 = self.get["mol1"]
        mol2 = self.get["mol2"]
        linked_mol = mol1.link(
            mol2, self["index1"], self["index2"], self.get["bond_order"]
        )
        if self.get("save_to_db", True):
            db.insert_derived_mol(
                linked_mol, update_duplicates=self.get("update_duplicates", False)
            )
        if self.get("save_mol_file", False):
            working_dir = os.getcwd()
            file_name = self.get(
                "filename.{}".format(self.get("fmt", "xyz")),
                "{}.{}".format(get_mol_formula(linked_mol), self.get("fmt", "xyz")),
            )
            linked_mol_file = os.path.join(working_dir, file_name)
            linked_mol.to(self.get("fmt", "xyz"), linked_mol_file)
        fw_spec["prev_calc_molecule"] = linked_mol


@explicit_serialize
class BreakMolecule(FiretaskBase):
    """
    credits: Samuel Blau
    """

    required_params = []
    optional_params = [
        "mol",
        "bonds",
        "ref_charge",
        "fragment_charges",
        "open_rings",
        "opt_steps",
        "working_dir",
        "db",
        "opt_gaussian_inputs",
        "freq_gaussian_inputs",
        "cart_coords",
        "oxidation_states",
        "calc_frags",
        "save_to_db",
        "save_to_file",
        "fmt",
        "update_duplicates",
        "additional_kwargs",
    ]

    @staticmethod
    def _define_charges(ref_charge, fragment_charges):
        # TODO: check charges on metal atoms so as to not violate valence rule
        # get a list of possible charges that each fragment can take
        possible_charges = []
        if ref_charge == 0:
            possible_charges.extend((-1, 0, 1))
        elif ref_charge > 0:
            for i in range(ref_charge + 1):
                possible_charges.append(ref_charge - i)
        else:
            for i in range(abs(ref_charge - 1)):
                possible_charges.append(ref_charge + i)
        # add additional charges to the list of possible charges
        if fragment_charges:
            fragment_charges += [ref_charge - charge for charge in fragment_charges]
            possible_charges.extend(
                charge for charge in fragment_charges if charge not in possible_charges
            )
        # find possible charge pairs upon breaking a bond; sum = ref_charge
        charge_pairs = [
            pair
            for pair in list(itertools.product(possible_charges, repeat=2))
            if sum(pair) == ref_charge
        ]

        charge_ind_map = {j: i for i, j in enumerate(possible_charges)}
        return possible_charges, charge_pairs, charge_ind_map

    @staticmethod
    def _find_unique_fragments(fragments_list):
        # get unique fragments from the list of all fragments
        # for each bond, store the new indexes of the unique fragments formed
        # upon splitting the molecule (can be used for analysis purposes later)
        unique_fragments = []
        fragments_indices = []
        for fragment in fragments_list:
            indices = []
            for i in fragment:
                flag = 0
                for ind, j in enumerate(unique_fragments):
                    if i.isomorphic_to(j):
                        flag = 1
                        indices.append(ind)
                        break
                if flag == 0:
                    indices.append(len(unique_fragments))
                    unique_fragments.append(i)
            fragments_indices.append(indices)
        return unique_fragments, fragments_indices

    @staticmethod
    def _find_unique_molecules(
        unique_fragments,
        fragment_charges,
        db,
        working_dir,
        save_to_db,
        update_duplicates,
        save_to_file,
        fmt,
        calc_frags,
    ):
        # create molecule objects from the unique fragments and set the charge
        # of each molecule
        molecules = [fragment.molecule for fragment in unique_fragments]
        unique_molecules = []
        if not calc_frags:
            # only saving if frags will not be calculated; otherwise saving
            # will be handled via the dynamically created fireworks; done to
            # avoid double saving
            if save_to_db:
                db = get_db(db) if db else get_db()
                for mol in molecules:
                    db.insert_molecule(mol, update_duplicates=update_duplicates)
            if save_to_file:
                for mol in molecules:
                    mol_formula = get_mol_formula(mol)
                    file = os.path.join(working_dir, f"{mol_formula}.{fmt}")
                    mol.to(fmt, file)

        for mol in molecules:
            for charge in fragment_charges:
                mol_copy = copy.deepcopy(mol)
                mol_copy.set_charge_and_spin(charge)
                unique_molecules.append(mol_copy)
        return unique_molecules

    @staticmethod
    def _find_molecule_indices(
        fragments_indices, fragment_charges, charge_ind_map, charge_pairs, offset=0
    ):
        num_charges = len(fragment_charges)
        molecule_indices = []
        for fragment in fragments_indices:
            split_possibilities = []
            for charge_pair in charge_pairs:
                split_possibility = [
                    offset + fragment[i] * num_charges + charge_ind_map[j]
                    for i, j in enumerate(charge_pair)
                ]
                split_possibilities.append(split_possibility)
            molecule_indices.append(split_possibilities)
        return molecule_indices

    def _cleanup_kwargs(self):
        additional_kwargs = self.get("additional_kwargs", {})
        kwargs = {
            i: j
            for i, j in additional_kwargs.items()
            if i
            not in self.required_params
            + self.optional_params
            + [
                "mol_operation_type",
                "dir_structure",
                "process_mol_func",
                "mol_name",
                "from_fw_spec",
                "skips",
            ]
        }
        return kwargs

    @staticmethod
    def _workflow(
        mol,
        gout_key,
        working_dir,
        db,
        opt_gaussian_inputs,
        freq_gaussian_inputs,
        cart_coords,
        oxidation_states,
        save_to_db,
        save_to_file,
        fmt,
        update_duplicates,
        **kwargs,
    ):

        from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS

        dir_structure = ["charge_{}".format(str(mol.charge))]
        mol_formula = get_mol_formula(mol)

        if len(mol) == 1:
            skips = ["opt"]
        else:
            skips = None

        job_name = "opt_freq"
        _, _, frag_fws = common_fw(
            mol_operation_type="get_from_mol",
            mol=mol,
            working_dir=working_dir,
            dir_structure=dir_structure,
            db=db,
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            process_mol_func=False,
            mol_name=mol_formula,
            from_fw_spec=False,
            skips=skips,
            gout_key=gout_key,
            save_to_db=save_to_db,
            save_to_file=save_to_file,
            fmt=fmt,
            update_duplicates=update_duplicates,
            **kwargs,
        )
        return Workflow(
            frag_fws,
            name="{}_{}".format(mol_formula, job_name),
            **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS},
        )

    def run_task(self, fw_spec):
        # get principle molecule from fw_spec or user input
        if fw_spec.get("prev_calc_molecule"):
            mol = fw_spec.get("prev_calc_molecule")
        elif self.get("molecule"):
            mol = self.get("molecule")
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )

        ref_charge = self.get("ref_charge", mol.charge)
        db = self.get("db")
        working_dir = self.get("working_dir", os.getcwd())
        calc_frags = self.get("calc_frags", False)
        save_to_file = self.get("save_to_file")
        save_to_db = self.get("save_to_db")
        fmt = self.get("fmt", "xyz")
        update_duplicates = self.get("update_duplicates", False)

        # break the bonds: either those specified by the user inputs or all
        # the bonds in the molecule; only supports breaking bonds or opening
        # ring in the principle molecule
        mol_graph = MoleculeGraph.with_local_env_strategy(
            mol, OpenBabelNN(), reorder=False, extend_structure=False
        )
        all_bonds = self.get("bonds", None)
        if not all_bonds:
            all_bonds = [
                tuple(sorted([idx1, idx2])) for idx1, idx2, _ in mol_graph.graph.edges
            ]

        fragments = []
        bonds = []
        ring_fragments = []
        ring_bonds = []

        # get fragments by splitting the molecule at each bond
        for idx, bond in enumerate(all_bonds):
            try:
                fragments.append(
                    mol_graph.split_molecule_subgraphs([bond], allow_reverse=True)
                )
                bonds.append(bond)
            except Exception as e:
                logger.info(e)
                if self.get("open_rings"):
                    logger.info(
                        "opening ring by breaking bond between "
                        "site {} and site {}".format(str(bond[0]), str(bond[1]))
                    )
                    ring_fragments.append(
                        [open_ring(mol_graph, [bond], self.get("opt_steps", 10000))]
                    )
                    ring_bonds.append(bond)
                else:
                    logger.info(
                        "encountered a ring bond; should set open_ring "
                        "to True to open the ring"
                    )

        unique_molecules = []
        molecule_indices = []
        if fragments:
            fragment_charges = self.get("fragment_charges", None)
            possible_charges, charge_pairs, charge_ind_map = self._define_charges(
                ref_charge, fragment_charges
            )
            unique_fragments, fragments_indices = self._find_unique_fragments(fragments)
            unique_molecules = self._find_unique_molecules(
                unique_fragments,
                possible_charges,
                db,
                working_dir,
                save_to_db,
                update_duplicates,
                save_to_file,
                fmt,
                calc_frags,
            )
            molecule_indices = self._find_molecule_indices(
                fragments_indices, possible_charges, charge_ind_map, charge_pairs
            )
        ring_unique_molecules = []
        ring_molecule_indices = []
        if ring_fragments:
            ring_unique_fragments, ring_fragments_indices = self._find_unique_fragments(
                ring_fragments
            )
            ring_unique_molecules = self._find_unique_molecules(
                ring_unique_fragments,
                [ref_charge],
                db,
                working_dir,
                save_to_db,
                update_duplicates,
                save_to_file,
                fmt,
                calc_frags,
            )
            charge_pairs = [(ref_charge,)]
            charge_ind_map = {ref_charge: 0}
            ring_molecule_indices = self._find_molecule_indices(
                ring_fragments_indices,
                [ref_charge],
                charge_ind_map,
                charge_pairs,
                len(unique_molecules),
            )
        all_molecules = unique_molecules + ring_unique_molecules

        update_spec = {
            "bonds": bonds + ring_bonds,
            "fragments": all_molecules,
            "molecule_indices": molecule_indices + ring_molecule_indices,
        }

        if calc_frags:
            wfs = []
            frag_keys = []
            opt_gaussian_inputs = self.get("opt_gaussian_inputs") or {}
            freq_gaussian_inputs = self.get("freq_gaussian_inputs") or {}
            cart_coords = self.get("cart_coords", True)
            oxidation_states = self.get("oxidation_states")
            additional_kwargs = self._cleanup_kwargs()
            for mol_ind, mol in enumerate(all_molecules):
                gout_key = "frag_{}".format(mol_ind)
                frag_wf = self._workflow(
                    mol,
                    gout_key,
                    working_dir,
                    db,
                    opt_gaussian_inputs,
                    freq_gaussian_inputs,
                    cart_coords,
                    oxidation_states,
                    save_to_db,
                    save_to_file,
                    fmt,
                    update_duplicates,
                    **additional_kwargs,
                )
                wfs.append(frag_wf)
                frag_keys.append(gout_key)
            update_spec["frag_keys"] = frag_keys
            return FWAction(update_spec=update_spec, detours=wfs, propagate=True)

        else:
            return FWAction(update_spec=update_spec)
