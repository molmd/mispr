from __future__ import annotations
from typing import Dict, Any
from . import defaults as D

def normalize_psi4_inputs(user: Dict[str, Any] | None) -> Dict[str, Any]:
    """Merge user params with robust defaults for Psi4 input writing."""
    u = dict(user or {})
    p = {
        "method": u.get("method", D.DEFAULT_METHOD),
        "basis": u.get("basis", D.DEFAULT_BASIS),
        "memory": u.get("memory", D.DEFAULT_MEMORY),

        # DEFAULT CHANGED: now opt+freq if user doesn't specify anything
        "job_type": _normalize_job_type(u.get("job_type", u.get("calculation", "opt+freq"))),

        # SCF/opt/grid controls
        "scf_type": u.get("scf_type", D.DEFAULT_SCF_TYPE),
        "e_convergence": u.get("e_convergence", D.DEFAULT_E_CONV),
        "d_convergence": u.get("d_convergence", D.DEFAULT_D_CONV),
        "g_convergence": u.get("g_convergence", D.DEFAULT_G_CONV),
        "dft_radial_points": u.get("dft_radial_points", D.DEFAULT_DFT_RADIAL_POINTS),
        "dft_spherical_points": u.get("dft_spherical_points", D.DEFAULT_DFT_SPHERICAL_POINTS),
        "maxiter": u.get("maxiter", D.DEFAULT_MAXITER),

        # reference: if user provided, keep; else we may auto-choose later
        "reference": u.get("reference", D.DEFAULT_REFERENCE),

        # optional extras pass-through: "options": {...} to add raw set-lines
        "options": u.get("options", None),

        # optional explicit charge/multiplicity overrides
        "charge": u.get("charge", None),
        "spin_multiplicity": u.get("spin_multiplicity", None),

        # optional threading: use -n on CLI typically; we allow set_num_threads() inline if desired
        "nthreads": u.get("nthreads", u.get("threads", u.get("n_cores", None))),

        # NEW DEFAULT: return the wavefunction from driver calls unless user disables
        "return_wfn": u.get("return_wfn", True),
    }
    return p

def _normalize_job_type(raw: str) -> str:
    r = (raw or "").lower()
    if r in ("opt", "optimize", "optimization"):
        return "opt"
    if r in ("freq", "frequency", "frequencies"):
        return "freq"
    if r in ("sp", "singlepoint", "single-point", "energy"):
        return "sp"
    if r in ("opt+freq", "freq+opt", "optfreq"):
        return "opt+freq"
    if r in ("gradient", "grad"):
        return "gradient"
    # default fallback
    return "sp"

def choose_reference(method: str, multiplicity: int) -> str:
    """Pick RHF/RKS for closed-shell, UHF/UKS for open-shell, based on method."""
    m = (method or "").lower()
    is_wf = m in ("hf", "scf", "mp2", "ccsd", "ccsd(t)", "qcisd", "mp3", "mp4")
    if multiplicity and int(multiplicity) > 1:
        return "uhf" if is_wf else "uks"
    # closed-shell: let Psi4 default to RHF/RKS unless user forced something
    return "rhf" if is_wf else "rks"

def build_set_block(params: Dict[str, Any]) -> str:
    """Return a ready-to-write 'set { ... }' block from normalized params."""
    # Assemble in a stable, readable order; omit None
    # Always include basis (this is the standard Psi4 pattern).
    kv = {
        "basis": params["basis"],
        "reference": params.get("reference", None),
        "scf_type": params.get("scf_type", None),
        "e_convergence": params.get("e_convergence", None),
        "d_convergence": params.get("d_convergence", None),
        "g_convergence": params.get("g_convergence", None),
        "maxiter": params.get("maxiter", None),
        "dft_radial_points": params.get("dft_radial_points", None),
        "dft_spherical_points": params.get("dft_spherical_points", None),
    }
    # allow raw `options` from user to append verbatim
    lines = ["set {"]
    for key in D.SETTABLE_KEYS_IN_ORDER:
        if key in kv and kv[key] is not None:
            lines.append(f"  {key} {kv[key]}")
    lines.append("}")
    # Append free-form options (simple 'set key value' per line)
    if isinstance(params.get("options"), dict):
        for k, v in params["options"].items():
            if k in kv:
                continue  # already written above
            v_str = "true" if v is True else "false" if v is False else str(v)
            lines.append(f"set {k} {v_str}")
    return "\n".join(lines)
