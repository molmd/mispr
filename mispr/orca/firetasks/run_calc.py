"""Define firetasks for running ORCA calculations.

ORCA follows the same external-process model as Gaussian (write an input file,
run the executable, parse the output file), but unlike Gaussian, there is no
pymatgen.io support for its file formats, so the input writing / output
parsing lives in mispr.orca.utilities. RunOrca combines the write/run/parse
steps into a single Firetask and stores a result dictionary with the same
schema as mispr.gaussian.utilities.gout.
process_run under fw_spec["gaussian_output"], so any downstream Firetask
written for Gaussian runs (BDEtoDB, ESPtoDB, BindingEnergytoDB, ...) can
consume ORCA results unmodified.

Error handling is currently a plain "did ORCA print its normal-termination
banner" check; wrapping the ORCA invocation in custodian (with ORCA-specific
error handlers, as mispr does for Gaussian via RunGaussianCustodian) is planned
follow-up work.
"""

import os
import json
import logging
import subprocess

from timeit import default_timer as timer

import numpy as np

from copy import deepcopy

from pymatgen.core.structure import Molecule, IMolecule

from fireworks.core.firework import FWAction, FiretaskBase
from fireworks.utilities.fw_utilities import explicit_serialize

from mispr.gaussian.defaults import JOB_TYPES
from mispr.gaussian.utilities.misc import recursive_signature_remove
from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.gaussian.utilities.db_utilities import get_db
from mispr.orca.utilities.inputs import get_orca_input_string
from mispr.orca.utilities.outputs import parse_orca_output

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"

logger = logging.getLogger(__name__)

DEFAULT_KEY = "gout_key"

# default ORCA job parameters used if not specified by the caller
DEFAULT_FUNCTIONAL = "b3lyp"
DEFAULT_BASIS_SET = "6-31G(d)"
DEFAULT_MEMORY_MB = 4000
DEFAULT_NUM_CORES = 1


GAS_CONSTANT = 8.314462618  # J/(mol K)
AVOGADRO = 6.02214076e23
BOLTZMANN = 1.380649e-23  # J/K
PLANCK = 6.62607015e-34  # J s
AMU_TO_KG = 1.66053906660e-27
HARTREE_TO_J = 4.3597447222071e-18
STANDARD_T = 298.15  # K
STANDARD_P = 101325.0  # Pa


def _json_default(o):
    """``json.dumps`` fallback for numpy scalars/arrays left over in a gout_dict; anything else is stringified."""
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _atomic_thermo_corrections(mass_amu, multiplicity, t=STANDARD_T, p=STANDARD_P):
    """
    Analytic ideal-gas thermochemistry corrections for a single free atom.

    A one-atom system has no vibrational degrees of freedom, so its frequency
    "calculation" -- which BDE workflows request for every fragment, including
    the lone H atoms produced by breaking X-H bonds -- reduces to the purely
    translational + electronic corrections, computed analytically here instead
    of invoking ORCA's Hessian machinery on a system it has nothing to do with.

    Returns:
        dict: "Zero-point correction", "Enthalpy", and "Gibbs Free Energy"
            corrections, in Hartree, using the same additive convention as
            mispr.gaussian (i.e. to be added to the raw electronic energy).

    """
    mass_kg = mass_amu * AMU_TO_KG
    q_trans = (2 * np.pi * mass_kg * BOLTZMANN * t / PLANCK**2) ** 1.5 * (
        BOLTZMANN * t / p
    )
    s_trans = GAS_CONSTANT * (np.log(q_trans) + 2.5)  # J/(mol K)
    s_elec = GAS_CONSTANT * np.log(multiplicity)  # J/(mol K)
    u_corr = 1.5 * GAS_CONSTANT * t  # J/mol; no ZPE, no rotational dof
    h_corr = u_corr + GAS_CONSTANT * t  # J/mol
    g_corr = h_corr - t * (s_trans + s_elec)  # J/mol

    def to_hartree(x_per_mol):
        return x_per_mol / AVOGADRO / HARTREE_TO_J

    return {
        "Zero-point correction": 0.0,
        "Enthalpy": to_hartree(h_corr),
        "Gibbs Free Energy": to_hartree(g_corr),
    }


def _resolve_orca_cmd(task):
    """Resolve the ORCA executable: the task's "orca_cmd" param, then the ORCA_CMD environment variable, then a bare "orca" (PATH lookup).

    Parallel runs (num_cores > 1) require the full absolute path -- an
    ORCA/OpenMPI requirement.
    """
    return task.get("orca_cmd") or os.environ.get("ORCA_CMD") or "orca"


def _run_orca(orca_cmd, input_path, output_path):
    """Run ORCA on input_path, streaming its stdout (where ORCA writes all its results) into output_path."""
    with open(output_path, "w") as f:
        # orca_cmd/input_path come from local job config, not untrusted input;
        # args are passed as a list (no shell) so there is no shell-injection
        # vector regardless.
        subprocess.run(
            [orca_cmd, input_path],
            stdout=f,
            stderr=subprocess.STDOUT,
            cwd=os.path.dirname(os.path.abspath(input_path)),
            shell=False,
        )


def _error_tail(output_path, n_lines=20):
    """Last n_lines of the output file, for error reporting on abnormal termination (ORCA prints its abort reason at the very end)."""
    try:
        with open(output_path) as f:
            return "\n".join(f.read().splitlines()[-n_lines:])
    except OSError:
        return f"(could not read {output_path})"


@explicit_serialize
class RunOrca(FiretaskBase):
    """Run an ORCA calculation (single point, geometry optimization, and/or frequency analysis) and store the result under fw_spec["gaussian_output"] using the same dictionary schema mispr uses for Gaussian runs.

    Args:
        molecule (Molecule, optional): pymatgen Molecule to run the calculation
            on; required unless ``prev_calc_key`` is provided.
        prev_calc_key (str, optional): Key of a previous run in
            fw_spec["gaussian_output"] whose output molecule should be used as
            the input structure (mirrors chaining an opt firework into a freq
            firework); ignored if ``molecule`` is provided.
        charge (int, optional): Molecular charge; defaults to the charge on
            ``molecule`` if it carries one, otherwise 0.
        multiplicity (int, optional): Spin multiplicity; defaults to the value
            on ``molecule`` if it carries one, otherwise 1. For multiplicity > 1
            ORCA switches to an unrestricted reference automatically.
        oxidation_states (dict, optional): Oxidation states used to derive the
            charge (e.g. {"Li": 1, "O": -2}).
        functional (str, optional): DFT functional (or "hf"/wavefunction method
            name) for the ``!`` keyword line; defaults to "b3lyp".
        basis_set (str, optional): Basis set, in ORCA's naming (e.g.
            "6-31G(d)", "def2-SVP"); defaults to "6-31G(d)".
        route_parameters (dict, optional): Gaussian-style route parameters used
            only to decide the job type; recognized keys (case insensitive) are
            "opt" and "freq" (both together produce an "Opt Freq" run), anything
            else falls back to a single point. Defaults to a single point.
        solvent (dict, optional): Implicit solvent options; supported key is
            "solvent" (an ORCA CPCM solvent name, e.g. "water"), producing a
            ``CPCM(<solvent>)`` keyword. If not provided, gas phase.
        memory (int, optional): Memory per core in MB (ORCA ``%maxcore``);
            defaults to 4000.
        num_cores (int, optional): Number of parallel processes (ORCA ``%pal``);
            defaults to 1. Values > 1 require ``orca_cmd`` to be an absolute
            path (ORCA/OpenMPI requirement).
        orca_cmd (str, optional): Path to the ORCA executable; falls back to the
            ORCA_CMD environment variable, then to "orca" on PATH.
        input_file (str, optional): Name of the input file to write; defaults to
            "mol.inp". The output file uses the same stem with ".out".
        db (str or dict, optional): Database credentials; path to db.json or a
            dict.
        save_to_db (bool, optional): Whether to insert the run into the runs
            collection.
        save_to_file (bool, optional): Whether to save the run to a json file.
        filename (str, optional): Name of the json file to save the run to;
            defaults to "run.json".
        gout_key (str, optional): Key to store this run under in
            fw_spec["gaussian_output"]; the run is always additionally stored
            under the default key "gout_key" (as with Gaussian runs).
        tag (str, optional): Tag stored in the db documents for easy retrieval.

    """

    required_params = []
    optional_params = [
        "molecule",
        "prev_calc_key",
        "charge",
        "multiplicity",
        "oxidation_states",
        "cart_coords",
        "functional",
        "basis_set",
        "route_parameters",
        "solvent",
        "memory",
        "num_cores",
        "orca_cmd",
        "input_file",
        "db",
        "save_to_db",
        "save_to_file",
        "filename",
        "gout_key",
        "tag",
    ]

    def _charge_from_oxidation_states(self, mol):
        """Calculate the charge of a molecule/cluster from the oxidation state of its individual elements (e.g. {"Li": 1, "O": -2}); mirrors mispr.gaussian.firetasks.write_inputs.WriteInput._update_charge, since this is plain pymatgen bookkeeping with no dependency on the QM engine."""
        mol_copy = deepcopy(mol)
        mol_copy.add_oxidation_state_by_element(self["oxidation_states"])
        mol_copy.set_charge_and_spin(super(IMolecule, mol_copy).charge)
        return int(mol_copy.charge)

    def _get_molecule(self, fw_spec):
        """Resolve the input molecule from (in priority order) the explicit "molecule" param, a previous run's optimized geometry via "prev_calc_key", or a linked molecule left in fw_spec by an earlier Firetask in the same Firework; raises KeyError if none is found."""
        mol = self.get("molecule")
        if mol is not None:
            if not isinstance(mol, Molecule):
                raise TypeError("molecule must be a pymatgen Molecule object")
            return mol

        prev_calc_key = self.get("prev_calc_key")
        if prev_calc_key:
            gout_dict = fw_spec.get("gaussian_output", {}).get(prev_calc_key)
            if not gout_dict:
                raise KeyError(
                    f"No previous run found under fw_spec['gaussian_output']"
                    f"['{prev_calc_key}']"
                )
            return Molecule.from_dict(gout_dict["output"]["output"]["molecule"])

        # set by mispr.gaussian.firetasks.geo_transformation.ProcessMoleculeInput
        # (e.g. after linking two previously-computed molecules into a complex,
        # ahead of this Firetask in the same Firework)
        if fw_spec.get("prev_calc_molecule"):
            return fw_spec["prev_calc_molecule"]

        raise KeyError(
            "No molecule present; provide 'molecule' or 'prev_calc_key', or "
            "check fw_spec"
        )

    def _resolve_charge_multiplicity(self, mol):
        """Resolve (charge, multiplicity) from the task params, falling back to mol's own values."""
        if self.get("oxidation_states"):
            charge = self._charge_from_oxidation_states(mol)
        else:
            charge = self.get("charge", getattr(mol, "charge", 0) or 0)

        if self.get("multiplicity") is not None:
            multiplicity = self.get("multiplicity")
        else:
            # the default multiplicity must be recomputed for whatever "charge"
            # ends up being (e.g. when it comes from oxidation_states, which can
            # differ from mol's own charge) -- mol's own spin_multiplicity
            # attribute reflects mol's original charge, not this one
            mol_for_mult = deepcopy(mol)
            mol_for_mult.set_charge_and_spin(charge)
            multiplicity = mol_for_mult.spin_multiplicity
        return charge, multiplicity

    def _build_job_spec(self, mol, charge):
        """Resolve the ORCA route keywords and opt/freq/PCM job flags from the task params."""
        functional = self.get("functional", DEFAULT_FUNCTIONAL)
        basis_set = self.get("basis_set", DEFAULT_BASIS_SET)
        route_parameters = self.get("route_parameters") or {}
        job_types = sorted(
            k.lower() for k in route_parameters if k.lower() in JOB_TYPES
        )
        if not job_types:
            job_types = ["sp"]

        solvent = self.get("solvent")
        is_pcm = bool(solvent)

        keywords = [functional, basis_set]
        run_opt = "opt" in job_types
        run_freq = "freq" in job_types

        # single-atom guards: an atom has no internal coordinates to optimize
        # and no vibrational modes to compute, so opt degenerates to a single
        # point and freq to analytic ideal-gas thermochemistry -- and a bare
        # nucleus (e.g. the H+ fragment BDE's charge-state enumeration
        # produces) has no electrons for ORCA's SCF to solve at all, so it is
        # skipped entirely with zero energy
        single_atom = mol.num_sites == 1
        bare_nucleus = single_atom and (mol.species[0].Z - charge) <= 0
        if single_atom:
            run_opt = False
        atomic_freq = run_freq and single_atom
        if atomic_freq:
            run_freq = False

        if run_opt:
            keywords.append("Opt")
        if run_freq:
            keywords.append("Freq")
        if is_pcm:
            keywords.append(f"CPCM({solvent.get('solvent', 'water')})")

        return {
            "functional": functional,
            "basis_set": basis_set,
            "route_parameters": route_parameters,
            "job_types": job_types,
            "keywords": keywords,
            "run_opt": run_opt,
            "run_freq": run_freq,
            "is_pcm": is_pcm,
            "bare_nucleus": bare_nucleus,
            "atomic_freq": atomic_freq,
        }

    def _extract_success_result(self, mol, charge, multiplicity, job_spec, parsed):
        """Pull energy/geometry/frequency/dipole/corrections out of a normally-terminated, converged ORCA run."""
        run_opt = job_spec["run_opt"]
        run_freq = job_spec["run_freq"]
        corrections = {}
        final_mol = mol
        frequencies = None

        energy = parsed["final_energy"]
        if run_opt and parsed["coords"]:
            final_mol = Molecule(parsed["species"], parsed["coords"])
            final_mol.set_charge_and_spin(charge, multiplicity)
        if run_freq:
            frequencies = parsed["frequencies"]
            # ORCA's thermochemistry section reports enthalpy/Gibbs as totals
            # (electronic energy + correction); mispr expects the correction
            # alone in all three cases, to be added back to final_energy
            # downstream (see BDEtoDB/IPEAtoDB)
            if parsed["zpe"] is not None:
                corrections["Zero-point correction"] = parsed["zpe"]
            if parsed["total_enthalpy"] is not None:
                corrections["Enthalpy"] = parsed["total_enthalpy"] - energy
            if parsed["gibbs_free_energy"] is not None:
                corrections["Gibbs Free Energy"] = (
                    parsed["gibbs_free_energy"] - energy
                )
        if job_spec["atomic_freq"]:
            corrections.update(
                _atomic_thermo_corrections(mol.species[0].atomic_mass, multiplicity)
            )
        dipole = parsed["dipole_moment"]

        return energy, final_mol, frequencies, corrections, dipole

    def _run_and_parse(self, mol, charge, multiplicity, job_spec, input_path, output_path):
        """Run ORCA (unless bare_nucleus) and collect the energy/geometry/frequency/error results."""
        run_opt = job_spec["run_opt"]
        corrections = {}
        has_completed = True
        error_msg = None
        energy = None
        final_mol = mol
        dipole = None
        frequencies = None
        orca_version = None

        if job_spec["bare_nucleus"]:
            energy = 0.0
            if job_spec["atomic_freq"]:
                corrections.update(
                    {
                        "Zero-point correction": 0.0,
                        "Enthalpy": 0.0,
                        "Gibbs Free Energy": 0.0,
                    }
                )
        else:
            with open(input_path, "w") as f:
                f.write(
                    get_orca_input_string(
                        mol,
                        charge,
                        multiplicity,
                        job_spec["keywords"],
                        memory_mb=self.get("memory", DEFAULT_MEMORY_MB),
                        num_cores=self.get("num_cores", DEFAULT_NUM_CORES),
                    )
                )
            _run_orca(_resolve_orca_cmd(self), input_path, output_path)
            parsed = parse_orca_output(output_path)
            orca_version = parsed["orca_version"]

            if not parsed["terminated_normally"]:
                has_completed = False
                error_msg = (
                    "ORCA did not terminate normally; end of output:\n"
                    + _error_tail(output_path)
                )
            elif run_opt and not parsed["opt_converged"]:
                has_completed = False
                error_msg = "ORCA geometry optimization did not converge"
            else:
                energy, final_mol, frequencies, corrections, dipole = (
                    self._extract_success_result(
                        mol, charge, multiplicity, job_spec, parsed
                    )
                )

        return {
            "energy": energy,
            "final_mol": final_mol,
            "dipole": dipole,
            "frequencies": frequencies,
            "corrections": corrections,
            "has_completed": has_completed,
            "error_msg": error_msg,
            "orca_version": orca_version,
        }

    def _build_gout_dict(self, mol, charge, multiplicity, job_spec, run_result, fw_spec, run_time):
        """Assemble the final gaussian_output-schema dict from the job spec and run result."""
        output_block = {
            "final_energy": run_result["energy"],
            "molecule": run_result["final_mol"].as_dict(),
        }
        if run_result["corrections"]:
            output_block["corrections"] = run_result["corrections"]
        if run_result["frequencies"]:
            output_block["frequencies"] = run_result["frequencies"]
        if run_result["dipole"]:
            output_block["dipole_moment"] = run_result["dipole"]
        if run_result["error_msg"]:
            output_block["error_message"] = run_result["error_msg"]

        orca_version = run_result["orca_version"]
        gout_dict = {
            "input": {
                "functional": job_spec["functional"],
                "basis_set": job_spec["basis_set"],
                "route_parameters": job_spec["route_parameters"],
                "charge": charge,
                "spin_multiplicity": multiplicity,
                "molecule": mol.as_dict(),
            },
            "output": {
                "output": output_block,
                "has_gaussian_completed": run_result["has_completed"],
                "is_pcm": job_spec["is_pcm"],
            },
            "functional": job_spec["functional"],
            "basis": job_spec["basis_set"],
            "phase": "solution" if job_spec["is_pcm"] else "gas",
            "type": ";".join(job_spec["job_types"]),
            **get_chem_schema(run_result["final_mol"]),
            "gauss_version": f"orca-{orca_version}" if orca_version else "orca",
        }
        gout_dict = {
            i: j
            for i, j in gout_dict.items()
            if i not in ["sites", "@module", "@class", "charge", "spin_multiplicity"]
        }
        if "tag" in fw_spec:
            gout_dict["tag"] = fw_spec["tag"]
        gout_dict["wall_time (s)"] = run_time
        gout_dict = json.loads(json.dumps(gout_dict, default=_json_default))
        return recursive_signature_remove(gout_dict)

    def _persist_result(self, working_dir, gout_dict):
        """Save gout_dict to db/file as requested and build the FWAction that stores it under fw_spec["gaussian_output"]."""
        run_list = {}
        db = self.get("db")
        if self.get("save_to_db"):
            runs_db = get_db(db)
            run_id = runs_db.insert_run(gout_dict)
            run_list["run_id_list"] = run_id
            logger.info("Saved parsed ORCA output to db")

        if self.get("save_to_file"):
            filename = self.get("filename", "run")
            file_path = os.path.join(working_dir, f"{filename}.json")
            with open(file_path, "w") as f:
                f.write(json.dumps(gout_dict, default=str))
            run_list["run_loc_list"] = file_path
            logger.info("Saved parsed ORCA output to json file")

        uid = self.get("gout_key")
        set_dict = {f"gaussian_output->{DEFAULT_KEY}": gout_dict}
        if uid:
            set_dict[f"gaussian_output->{uid}"] = gout_dict
        mod_dict = {"_set": set_dict}
        if run_list:
            mod_dict.update({"_push": run_list})
        return FWAction(mod_spec=mod_dict, propagate=True)

    def run_task(self, fw_spec):
        """Write the ORCA input file, run ORCA, parse its output, and store the result under fw_spec["gaussian_output"] for downstream Firetasks (e.g. ``mispr.gaussian.firetasks.parse_outputs.ESPtoDB``) to consume."""
        working_dir = os.getcwd()
        mol = self._get_molecule(fw_spec)
        charge, multiplicity = self._resolve_charge_multiplicity(mol)

        if self.get("cart_coords") is False:
            raise NotImplementedError(
                "The ORCA backend only writes cartesian-coordinate inputs; "
                "z-matrix input (cart_coords=False) is not supported"
            )

        job_spec = self._build_job_spec(mol, charge)

        input_file = self.get("input_file", "mol.inp")
        input_path = os.path.join(working_dir, input_file)
        output_path = os.path.splitext(input_path)[0] + ".out"

        st = timer()
        run_result = self._run_and_parse(
            mol, charge, multiplicity, job_spec, input_path, output_path
        )
        run_time = timer() - st
        fw_spec["run_time"] = run_time

        gout_dict = self._build_gout_dict(
            mol, charge, multiplicity, job_spec, run_result, fw_spec, run_time
        )

        if not run_result["has_completed"]:
            raise ValueError(f"ORCA did not complete normally: {run_result['error_msg']}")

        return self._persist_result(working_dir, gout_dict)


@explicit_serialize
class ESP(FiretaskBase):
    """Compute ESP-fitted atomic partial charges for a molecule via ORCA's CHELPG scheme (grid-based ESP fitting -- ORCA's counterpart to the Merz-Singh-Kollman fit the Gaussian backend uses; both fit atomic charges to the molecular electrostatic potential, differing in grid construction/restraints).

    Always runs a single-point calculation -- callers are expected to have
    already optimized the molecule (e.g. via a preceding ``RunOrca``
    Firework) and pass it in through ``prev_calc_key``.

    Args:
        molecule (Molecule, optional): pymatgen Molecule to run the ESP
            calculation on directly; mutually exclusive with ``prev_calc_key``
            (one of the two must be provided).
        prev_calc_key (str, optional): Key into fw_spec["gaussian_output"] whose
            optimized geometry should be used as the input structure.
        charge (int, optional): Charge on the molecule; defaults to the
            molecule's own charge (or 0).
        multiplicity (int, optional): Spin multiplicity; defaults to the
            molecule's own value (or 1).
        method_esp (str, optional): Method for the ESP single-point calculation;
            defaults to "hf".
        basis_esp (str, optional): Basis set for the ESP calculation, in ORCA's
            naming; defaults to "6-31G*".
        memory (int, optional): Memory per core in MB; defaults to 4000.
        num_cores (int, optional): Number of parallel processes; defaults to 1.
        orca_cmd (str, optional): Path to the ORCA executable; falls back to the
            ORCA_CMD environment variable, then to "orca" on PATH.
        input_file (str, optional): Name of the input file to write; defaults to
            "esp.inp".
        db (str or dict, optional): Database credentials; path to db.json or a
            dict.
        save_to_db (bool, optional): Whether to insert the run into the runs
            collection.
        save_to_file (bool, optional): Whether to save the run to a json file.
        filename (str, optional): Name of the json file to save the run to.
        gout_key (str, optional): Key to store this run under in
            fw_spec["gaussian_output"]; the run is always additionally stored
            under the default key "gout_key" (as with Gaussian runs).
        tag (str, optional): Tag stored in the db documents for easy retrieval.

    """

    required_params = []
    optional_params = [
        "molecule",
        "prev_calc_key",
        "charge",
        "multiplicity",
        "method_esp",
        "basis_esp",
        "memory",
        "num_cores",
        "orca_cmd",
        "input_file",
        "db",
        "save_to_db",
        "save_to_file",
        "filename",
        "gout_key",
        "tag",
    ]

    def run_task(self, fw_spec):
        """Run the CHELPG ESP single-point calculation, then store the result under fw_spec["gaussian_output"] for a downstream ``mispr.gaussian.firetasks.parse_outputs.ESPtoDB`` to consume."""
        working_dir = os.getcwd()
        mol = self.get("molecule")
        if mol is not None:
            pass
        elif self.get("prev_calc_key"):
            prev_calc_key = self.get("prev_calc_key")
            mol = Molecule.from_dict(
                fw_spec.get("gaussian_output", {}).get(prev_calc_key)["output"][
                    "output"
                ]["molecule"]
            )
        else:
            raise KeyError("Don't have 'molecule' and 'prev_calc_key' ")

        charge = self.get("charge", getattr(mol, "charge", 0) or 0)
        multiplicity = self.get(
            "multiplicity", getattr(mol, "spin_multiplicity", 1) or 1
        )

        method_esp = self.get("method_esp", "hf")
        basis_esp = self.get("basis_esp", "6-31G*")
        keywords = [method_esp, basis_esp, "CHELPG"]

        input_file = self.get("input_file", "esp.inp")
        input_path = os.path.join(working_dir, input_file)
        output_path = os.path.splitext(input_path)[0] + ".out"

        has_completed = True
        error_msg = None
        energy = None
        esp_charges = None
        orca_version = None

        with open(input_path, "w") as f:
            f.write(
                get_orca_input_string(
                    mol,
                    charge,
                    multiplicity,
                    keywords,
                    memory_mb=self.get("memory", DEFAULT_MEMORY_MB),
                    num_cores=self.get("num_cores", DEFAULT_NUM_CORES),
                )
            )
        _run_orca(_resolve_orca_cmd(self), input_path, output_path)
        parsed = parse_orca_output(output_path)
        orca_version = parsed["orca_version"]

        if not parsed["terminated_normally"]:
            has_completed = False
            error_msg = (
                "ORCA did not terminate normally; end of output:\n"
                + _error_tail(output_path)
            )
        elif parsed["chelpg_charges"] is None:
            has_completed = False
            error_msg = "no CHELPG charges found in ORCA output"
        else:
            energy = parsed["final_energy"]
            esp_charges = parsed["chelpg_charges"]

        output_block = {
            "final_energy": energy,
            "molecule": mol.as_dict(),
        }
        if esp_charges is not None:
            output_block["ESP_charges"] = esp_charges
        if error_msg:
            output_block["error_message"] = error_msg

        gout_dict = {
            "input": {
                "functional": method_esp,
                "basis_set": basis_esp,
                "charge": charge,
                "spin_multiplicity": multiplicity,
                "molecule": mol.as_dict(),
            },
            "output": {
                "output": output_block,
                "has_gaussian_completed": has_completed,
                "is_pcm": False,
            },
            "functional": method_esp,
            "basis": basis_esp,
            "phase": "gas",
            "type": "esp",
            **get_chem_schema(mol),
            "gauss_version": f"orca-{orca_version}" if orca_version else "orca",
        }
        gout_dict = {
            i: j
            for i, j in gout_dict.items()
            if i not in ["sites", "@module", "@class", "charge", "spin_multiplicity"]
        }
        if "tag" in fw_spec:
            gout_dict["tag"] = fw_spec["tag"]
        gout_dict = json.loads(json.dumps(gout_dict, default=_json_default))
        gout_dict = recursive_signature_remove(gout_dict)

        if not has_completed:
            raise ValueError(
                f"ORCA ESP calculation did not complete normally: {error_msg}"
            )

        run_list = {}
        db = self.get("db")
        if self.get("save_to_db"):
            runs_db = get_db(db)
            run_id = runs_db.insert_run(gout_dict)
            run_list["run_id_list"] = run_id
            logger.info("Saved parsed ORCA ESP output to db")

        if self.get("save_to_file"):
            filename = self.get("filename", "run")
            file_path = os.path.join(working_dir, f"{filename}.json")
            with open(file_path, "w") as f:
                f.write(json.dumps(gout_dict, default=str))
            run_list["run_loc_list"] = file_path
            logger.info("Saved parsed ORCA ESP output to json file")

        uid = self.get("gout_key")
        set_dict = {f"gaussian_output->{DEFAULT_KEY}": gout_dict}
        if uid:
            set_dict[f"gaussian_output->{uid}"] = gout_dict
        mod_dict = {"_set": set_dict}
        if run_list:
            mod_dict.update({"_push": run_list})
        return FWAction(mod_spec=mod_dict, propagate=True)
