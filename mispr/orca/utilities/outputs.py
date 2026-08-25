"""Parse ORCA output files.

ORCA reports everything through a single free-text output stream (what the
``orca`` binary prints to stdout), so -- in the role pymatgen.io.gaussian's
GaussianOutput plays for Gaussian -- this module extracts the quantities mispr's
workflows need (final energy, optimized geometry, vibrational frequencies,
thermochemistry, CHELPG charges, dipole moment) from that text.

All energies are returned in Hartree (ORCA's native output unit), coordinates in
Angstrom, frequencies in cm**-1, dipole components in atomic units -- matching
the conventions of the Gaussian backend so downstream analysis firetasks
(BDEtoDB, ESPtoDB, BindingEnergytoDB) can consume the results unmodified.
"""

import re

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"

# ORCA prints this banner (only) when a run finishes without error; its absence
# is the reliable failure signal -- the process exit code is not trustworthy
# for this (ORCA can exit 0 after an aborted calculation)
NORMAL_TERMINATION = "ORCA TERMINATED NORMALLY"
OPT_CONVERGED = "THE OPTIMIZATION HAS CONVERGED"

_FINAL_ENERGY_RE = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
_VERSION_RE = re.compile(r"Program Version\s+(\S+)")
_COORD_HEADER = "CARTESIAN COORDINATES (ANGSTROEM)"
_COORD_LINE_RE = re.compile(
    r"^\s*([A-Za-z]{1,2})\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$"
)
_FREQ_LINE_RE = re.compile(r"^\s*\d+:\s+(-?\d+\.\d+)\s+cm\*\*-1")
_ZPE_RE = re.compile(r"Zero point energy\s+\.\.\.\s+(-?\d+\.\d+)\s+Eh")
_ENTHALPY_RE = re.compile(r"Total Enthalpy\s+\.\.\.\s+(-?\d+\.\d+)\s+Eh")
_GIBBS_RE = re.compile(r"Final Gibbs free energy\s+\.\.\.\s+(-?\d+\.\d+)\s+Eh")
_DIPOLE_RE = re.compile(
    r"Total Dipole Moment\s+:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)"
)
_CHELPG_HEADER = "CHELPG Charges"
_CHELPG_LINE_RE = re.compile(r"^\s*\d+\s+[A-Za-z]{1,2}\s*:\s+(-?\d+\.\d+)")


def _parse_last_coord_block(lines):
    """Return (species, coords) from the last "CARTESIAN COORDINATES (ANGSTROEM)" block in the output -- for an optimization that is the converged geometry, since ORCA reprints the block after every step."""
    last_start = None
    for i, line in enumerate(lines):
        if _COORD_HEADER in line:
            last_start = i
    if last_start is None:
        return None

    species, coords = [], []
    # the header is followed by one dashed separator line, then the atom lines
    for line in lines[last_start + 2 :]:
        match = _COORD_LINE_RE.match(line)
        if not match:
            break
        species.append(match.group(1))
        coords.append([float(match.group(2)), float(match.group(3)), float(match.group(4))])
    if not species:
        return None
    return species, coords


def _parse_frequencies(lines):
    """Return the frequency list (cm**-1) from the last "VIBRATIONAL FREQUENCIES" section, including the leading zero entries ORCA prints for the translational/rotational modes (callers filter those out; imaginary modes appear as negative values)."""
    last_start = None
    for i, line in enumerate(lines):
        if "VIBRATIONAL FREQUENCIES" in line:
            last_start = i
    if last_start is None:
        return None

    freqs = []
    started = False
    for line in lines[last_start + 1 :]:
        match = _FREQ_LINE_RE.match(line)
        if match:
            started = True
            freqs.append(float(match.group(1)))
        elif started:
            break
    return freqs or None


def _parse_chelpg_charges(lines):
    """Return the per-atom charge list from the last "CHELPG Charges" block."""
    last_start = None
    for i, line in enumerate(lines):
        if _CHELPG_HEADER in line:
            last_start = i
    if last_start is None:
        return None

    charges = []
    started = False
    for line in lines[last_start + 1 :]:
        match = _CHELPG_LINE_RE.match(line)
        if match:
            started = True
            charges.append(float(match.group(1)))
        elif started:
            break
    return charges or None


def parse_orca_output(file_path):
    """
    Parse an ORCA output file into a plain dict of the quantities mispr needs.

    Every key is always present; a quantity the run did not produce (e.g.
    frequencies in a plain single point) is None. Where ORCA prints a quantity
    multiple times over the course of a run (energies, geometries, dipoles --
    once per optimization step), the last occurrence is returned, since that is
    the one belonging to the final geometry.

    Args:
        file_path (str): Path to the ORCA output file.

    Returns:
        dict: With keys:

            * "terminated_normally" (bool): Whether the run printed ORCA's
              normal-termination banner (the reliable success signal; the
              process exit code is not).
            * "orca_version" (str or None): ORCA version from the file header.
            * "final_energy" (float or None): Last "FINAL SINGLE POINT ENERGY",
              in Hartree.
            * "opt_converged" (bool): Whether a geometry optimization reported
              convergence (False for non-opt runs).
            * "species" (list of str or None) / "coords" (list of [x, y, z] or
              None): Last printed cartesian geometry, in Angstrom.
            * "frequencies" (list of float or None): Vibrational frequencies in
              cm**-1, zero translational/rotational entries removed; imaginary
              modes appear as negative values.
            * "zpe" (float or None): Zero point energy correction, in Hartree.
            * "total_enthalpy" (float or None): Total enthalpy (electronic
              energy + correction), in Hartree.
            * "gibbs_free_energy" (float or None): Final Gibbs free energy
              (electronic energy + correction), in Hartree.
            * "chelpg_charges" (list of float or None): CHELPG per-atom charges.
            * "dipole_moment" (list of float or None): Total dipole moment
              [x, y, z], in atomic units.

    """
    with open(file_path) as f:
        text = f.read()
    lines = text.splitlines()

    result = {
        "terminated_normally": NORMAL_TERMINATION in text,
        "orca_version": None,
        "final_energy": None,
        "opt_converged": OPT_CONVERGED in text,
        "species": None,
        "coords": None,
        "frequencies": None,
        "zpe": None,
        "total_enthalpy": None,
        "gibbs_free_energy": None,
        "chelpg_charges": None,
        "dipole_moment": None,
    }

    version_match = _VERSION_RE.search(text)
    if version_match:
        result["orca_version"] = version_match.group(1)

    energy_matches = _FINAL_ENERGY_RE.findall(text)
    if energy_matches:
        result["final_energy"] = float(energy_matches[-1])

    coord_block = _parse_last_coord_block(lines)
    if coord_block:
        result["species"], result["coords"] = coord_block

    freqs = _parse_frequencies(lines)
    if freqs is not None:
        # ORCA lists all 3N modes, printing exact 0.00 for the 5/6
        # translational + rotational ones; only the true vibrations are of
        # interest downstream
        result["frequencies"] = [f for f in freqs if f != 0.0]

    for key, regex in [
        ("zpe", _ZPE_RE),
        ("total_enthalpy", _ENTHALPY_RE),
        ("gibbs_free_energy", _GIBBS_RE),
    ]:
        matches = regex.findall(text)
        if matches:
            result[key] = float(matches[-1])

    dipole_matches = _DIPOLE_RE.findall(text)
    if dipole_matches:
        result["dipole_moment"] = [float(x) for x in dipole_matches[-1]]

    result["chelpg_charges"] = _parse_chelpg_charges(lines)

    return result
