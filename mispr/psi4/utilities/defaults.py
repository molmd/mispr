"""
Default settings for writing Psi4 input files (DFT-first, user-overridable).
These are *input-writing defaults* only; Psi4 itself will still decide runtime
details unless you set them here explicitly.
"""

# Method/basis & resources
DEFAULT_METHOD = "b3lyp"         # user can override (supports HF/MP2/CCSD/etc.)
DEFAULT_BASIS = "6-31g*"         # common starter basis; override freely
DEFAULT_MEMORY = "2 GB"          # appears as: memory 2 GB
DEFAULT_SCF_TYPE = "df"          # density-fitting where available

# Tight SCF & geometry convergence (Gaussian-like)
DEFAULT_E_CONV = 1e-8            # energy convergence
DEFAULT_D_CONV = 1e-8            # density convergence
DEFAULT_G_CONV = "gau_tight"     # optking geometry convergence preset

# DFT grid (Psi4 manual commonly recommends 99/590 for robust DFT)
DEFAULT_DFT_RADIAL_POINTS = 99
DEFAULT_DFT_SPHERICAL_POINTS = 590

# Reasonable other knobs (left None => omitted)
DEFAULT_MAXITER = None           # e.g., 200 to force in set{}
DEFAULT_REFERENCE = None         # auto-chosen from multiplicity unless user sets

# File names
DEFAULT_INPUT_FILENAME = "psi4_input.dat"
DEFAULT_OUTPUT_FILENAME = "psi4_output.out"

# Keys we’ll write inside `set {}` when present
# (Order chosen to be readable; omit if value is None)
SETTABLE_KEYS_IN_ORDER = (
    "basis",
    "reference",
    "scf_type",
    "e_convergence",
    "d_convergence",
    "g_convergence",
    "maxiter",
    "dft_radial_points",
    "dft_spherical_points",
    # add more here if you later want them emitted by default
)
