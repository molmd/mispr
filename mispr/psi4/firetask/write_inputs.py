import os
from copy import deepcopy

# --- pymatgen (match Gaussian writer style) ---
from pymatgen.core import Molecule
try:
    from pymatgen.core.structure import IMolecule  # older/newer versions differ
except Exception:
    IMolecule = Molecule  # safe fallback

# Optional graph-based isomorphism check
try:
    from pymatgen.analysis.local_env import OpenBabelNN
    from pymatgen.analysis.graphs import MoleculeGraph
    _HAS_OBABEL = True
except Exception:
    _HAS_OBABEL = False

# --- FireWorks ---
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize

# --- local utils ---
from mispr.psi4.utilities import defaults as D
from mispr.psi4.utilities.inputs import (
    normalize_psi4_inputs,
    choose_reference,
    build_set_block,
)

@explicit_serialize
class WritePsi4Input(FiretaskBase):
    """
    Write a Psi4 input (.dat) file.

    Optional params:
        molecule (pymatgen Molecule/IMolecule or as_dict):
            If absent, will look for fw_spec['prev_calc_molecule'] (dict or Molecule).
        psi4_input_params (dict):
            method, basis, memory, job_type ('sp'|'opt'|'freq'|'opt+freq'|'gradient'),
            scf_type, e_convergence, d_convergence, g_convergence,
            dft_radial_points, dft_spherical_points,
            reference (else auto-chosen from multiplicity),
            charge, spin_multiplicity (override Molecule’s),
            options: { arbitrary_key: value }  # appended as 'set key value'
            nthreads: int  # emits set_num_threads(n)
            return_wfn: bool  # default True
        input_file (str): filename to write (default: psi4_input.dat)
        output_file (str): recorded to fw_spec for later runners (default: psi4_output.out)
        oxidation_states (dict): e.g. {'Na': +1, 'Cl': -1} to infer total charge
        cart_coords (bool): reserved (Cartesian always written for now)
        set_reference_from_mult (bool): default True (auto-pick RHF/RKS vs UHF/UKS)
    """
    required_params = []
    optional_params = [
        "molecule",
        "psi4_input_params",
        "input_file",
        "output_file",
        "oxidation_states",
        "cart_coords",
        "set_reference_from_mult",
    ]

    # ------------------------------ helpers ------------------------------

    def _coerce_mol(self, obj):
        """Accept Molecule/IMolecule/dict and return a Molecule."""
        if isinstance(obj, Molecule):
            return obj
        if isinstance(obj, IMolecule):
            return Molecule(obj.species, obj.cart_coords, charge=getattr(obj, "charge", 0),
                            spin_multiplicity=getattr(obj, "spin_multiplicity", 1))
        if isinstance(obj, dict):
            try:
                return Molecule.from_dict(obj)
            except Exception:
                pass
        raise ValueError("WritePsi4Input needs a pymatgen Molecule/IMolecule or as_dict() form.")

    def _maybe_apply_oxidation_charge(self, mol: Molecule):
        """If oxidation_states given, set formal-charge-based charge on Molecule."""
        ox = self.get("oxidation_states")
        if not ox:
            return
        tmp = deepcopy(mol)
        tmp.add_oxidation_state_by_element(ox)
        total = 0
        for site in tmp.sites:
            sp = site.specie
            if hasattr(sp, "oxi_state"):
                total += sp.oxi_state
        mol.set_charge_and_spin(int(total), mol.spin_multiplicity)

    def _select_geometry(self, fw_spec) -> Molecule:
        """Choose geometry from params vs fw_spec, prefer prev_calc_molecule if isomorphic."""
        mol_param = self.get("molecule", None)
        mol_spec = fw_spec.get("prev_calc_molecule", None)

        if mol_param is None and mol_spec is None:
            raise ValueError("No molecule provided: pass 'molecule' or provide fw_spec['prev_calc_molecule'].")

        if mol_param and mol_spec:
            m1 = self._coerce_mol(mol_param)
            m2 = self._coerce_mol(mol_spec)
            if _HAS_OBABEL:
                try:
                    g1 = MoleculeGraph.with_local_env_strategy(m1, OpenBabelNN(), reorder=False, extend_structure=False)
                    g2 = MoleculeGraph.with_local_env_strategy(m2, OpenBabelNN(), reorder=False, extend_structure=False)
                    if g1.isomorphic_to(g2):
                        return m2
                except Exception:
                    pass
            # If graphing fails or non-isomorphic, fall back to composition check
            if m1.composition.almost_equals(m2.composition):
                return m2
            return m1

        return self._coerce_mol(mol_param or mol_spec)

    # ------------------------------ main ------------------------------

    def run_task(self, fw_spec):
        # 1) pick geometry
        mol = self._select_geometry(fw_spec)

        # 2) optionally set charge from oxidation states
        self._maybe_apply_oxidation_charge(mol)

        # 3) merge params with defaults & allow explicit charge/mult override
        params = normalize_psi4_inputs(self.get("psi4_input_params"))
        if params.get("charge") is not None or params.get("spin_multiplicity") is not None:
            mol.set_charge_and_spin(
                int(params.get("charge", mol.charge or 0)),
                int(params.get("spin_multiplicity", mol.spin_multiplicity or 1)),
            )

        # 4) determine reference if requested
        if self.get("set_reference_from_mult", True) and not params.get("reference"):
            params["reference"] = choose_reference(params["method"], mol.spin_multiplicity)

        # 5) file names
        iname = self.get("input_file", D.DEFAULT_INPUT_FILENAME)
        oname = self.get("output_file", D.DEFAULT_OUTPUT_FILENAME)
        ipath = os.path.join(os.getcwd(), iname)

        # 6) write the Psi4 input file plainly (no post-processing gymnastics)
        with open(ipath, "w") as f:
            # memory
            mem = params.get("memory")
            if mem:
                f.write(f"memory {mem}\n\n")

            # optional set_num_threads
            if params.get("nthreads"):
                try:
                    n = int(params["nthreads"])
                    if n > 0:
                        f.write(f"set_num_threads({n})\n\n")
                except Exception:
                    pass

            # molecule block (Cartesian)
            q = int(getattr(mol, "charge", 0) or 0)
            mult = int(getattr(mol, "spin_multiplicity", 1) or 1)
            f.write("molecule {\n")
            f.write(f"  {q} {mult}\n")
            for sp, (x, y, z) in zip(mol.species, mol.cart_coords):
                f.write(f"  {sp.symbol:<2} {x:.8f} {y:.8f} {z:.8f}\n")
            f.write("}\n\n")

            # set-block with basis + defaults
            f.write(build_set_block({**params, "basis": params["basis"]}) + "\n\n")

            # driver calls (DEFAULT includes return_wfn=True now)
            method = params["method"].lower()
            job = params["job_type"]
            ret = params.get("return_wfn", True)

            def emit(fun: str):
                if ret:
                    f.write(f"{fun}('{method}', return_wfn=True)\n")
                else:
                    f.write(f"{fun}('{method}')\n")

            if job == "opt":
                emit("optimize")
            elif job == "freq":
                emit("frequency")
            elif job == "opt+freq":
                emit("optimize")
                emit("frequency")
            elif job == "gradient":
                emit("gradient")
            else:  # "sp"
                emit("energy")

        # 7) update spec (so later run/parse tasks know names)
        uspec = dict(fw_spec)
        uspec["psi4_input_file"] = iname
        uspec["psi4_output_file"] = oname
        uspec["psi4_method"] = params["method"]
        uspec["psi4_basis"] = params["basis"]
        uspec["psi4_reference"] = params.get("reference")
        return FWAction(update_spec=uspec)
