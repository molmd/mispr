"""Build ORCA input files from pymatgen Molecule objects.

ORCA runs as an external process reading a text input file, so -- exactly as
pymatgen.io.gaussian's GaussianInput does for Gaussian -- the pymatgen
Molecule has to be serialized into ORCA's input format first:

    ! B3LYP 6-31G(d) Opt Freq        <- "keyword line": method, basis, job types
    %maxcore 4000                    <- memory per core, in MB
    %pal nprocs 4 end                <- parallelism (omitted when serial)

    * xyz 0 1                        <- charge, spin multiplicity
    O    0.0000000000    0.0000000000    0.1173000000
    H    0.0000000000    0.7572000000   -0.4692000000
    H    0.0000000000   -0.7572000000   -0.4692000000
    *
"""

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"


def get_orca_input_string(
    mol,
    charge,
    multiplicity,
    keywords,
    memory_mb=4000,
    num_cores=1,
    extra_blocks=None,
):
    """
    Serialize a pymatgen Molecule into a complete ORCA input file string.

    Args:
        mol (Molecule): pymatgen Molecule to write (cartesian coordinates only --
            ORCA's z-matrix input uses a different, Gaussian-incompatible block
            format that mispr does not need).
        charge (int): Molecular charge, written into the ``* xyz`` line.
        multiplicity (int): Spin multiplicity, written into the ``* xyz`` line.
            For multiplicity > 1 ORCA switches to an unrestricted (UHF/UKS)
            reference automatically, so no extra keyword is needed for radicals.
        keywords (list of str): Entries for the ``!`` keyword line, in order --
            typically [functional, basis_set, job type keywords...]. Names must
            be ones ORCA recognizes (e.g. "B3LYP", "6-31G(d)", "Opt", "Freq",
            "CHELPG", "CPCM(water)").
        memory_mb (int, optional): Memory per core in MB (ORCA's ``%maxcore``);
            defaults to 4000.
        num_cores (int, optional): Number of parallel processes (``%pal``);
            defaults to 1 (serial), in which case the ``%pal`` block is omitted
            entirely. Note that running ORCA in parallel requires invoking the
            ``orca`` binary by its full absolute path (an ORCA/OpenMPI
            requirement, not a mispr one).
        extra_blocks (list of str, optional): Raw ``%...`` input blocks appended
            verbatim after the standard ones, for options with no keyword-line
            shorthand.

    Returns:
        str: The full ORCA input file content.

    """
    lines = ["! " + " ".join(str(k) for k in keywords)]
    lines.append(f"%maxcore {int(memory_mb)}")
    if num_cores and int(num_cores) > 1:
        lines.append(f"%pal nprocs {int(num_cores)} end")
    for block in extra_blocks or []:
        lines.append(block)
    lines.append("")
    lines.append(f"* xyz {int(charge)} {int(multiplicity)}")
    for site in mol:
        lines.append(
            f"{site.specie.symbol:<3} {site.x:>15.10f} {site.y:>15.10f} "
            f"{site.z:>15.10f}"
        )
    lines.append("*")
    lines.append("")
    return "\n".join(lines)
