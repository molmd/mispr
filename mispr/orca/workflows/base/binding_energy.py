"""Define the binding energy workflow, using ORCA instead of Gaussian.

Mirrors ``mispr.gaussian.workflows.base.binding_energy.get_binding_energies``.

* **Fireworks 1-4**: Optimize + frequency for each of the two molecules.
* **Firework 5**: Combine the two optimized molecules at ``index`` and optimize
  the resulting complex (``LinkedMolOrcaFW`` -- linking and optimizing must
  happen in the same Firework, since the linked molecule is only passed via
  ``fw_spec["prev_calc_molecule"]``, which does not survive a Firework
  boundary).
* **Firework 6**: Frequency calculation on the optimized complex.
* **Firework 7**: Compute the binding energy and save it (``BindingEnergytoDB``,
  reused unmodified from the Gaussian firetasks -- pure energy bookkeeping, not
  tied to any QM engine).

No counterpoise (BSSE) correction option is offered yet: it needs ghost-atom
support in the ORCA input writer (ORCA marks ghost atoms with a ":" suffix on
the element symbol), which is left as follow-up work.
"""

import os

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.firetasks.parse_outputs import BindingEnergytoDB
from mispr.orca.fireworks.core import OrcaFW, LinkedMolOrcaFW
from mispr.orca.workflows.base.core import common_fw, WORKFLOW_KWARGS

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"


def _to_orca_solvent(solvent_gaussian_inputs):
    """Translate the Gaussian-style solvent string (e.g. "(Solvent=Water)") into the dict RunOrca expects (e.g. {"solvent": "water"}); returns None for gas phase."""
    if not solvent_gaussian_inputs:
        return None
    solvent_inputs = [
        i.lower() for i in solvent_gaussian_inputs.strip("()").split(",")
    ]
    solvent_name = next(
        (s.split("=")[1] for s in solvent_inputs if "solvent" in s), "water"
    )
    return {"solvent": solvent_name}


def get_binding_energies(
    mol_operation_type,
    mol,
    index,
    bond_order=1,
    db=None,
    name="binding_energy_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    orca_cmd=None,
    num_cores=None,
    memory=None,
    **kwargs,
):
    """Define a workflow for calculating the binding energy between two molecules using ORCA: optimize + frequency each of the two molecules separately, link them together at the given atom indices and optimize the resulting complex, run its frequency, then compute the binding energy.

    Mirrors ``mispr.gaussian.workflows.base.binding_energy.get_binding_energies``,
    which has the same overall structure and arguments.

    ``opt_gaussian_inputs``/``freq_gaussian_inputs`` use the same keys as the
    Gaussian workflow (e.g. {"functional": "B3LYP", "basis_set": "6-31G(d)",
    "route_parameters": {"Opt": None}}); Gaussian-only keys are accepted but
    ignored, and basis set names must be ones ORCA recognizes.
    ``solvent_gaussian_inputs``/``solvent_properties`` keep their Gaussian-style
    format for ``BindingEnergytoDB``'s db metadata bookkeeping; they are
    translated internally into the CPCM option RunOrca expects to actually
    drive the implicit-solvent calculation.

    Args:
        mol_operation_type (list): Two-element list of molecule operation
            types, one per molecule; see ``process_mol`` in
            ``mispr/gaussian/utilities/mol.py`` for supported operations.
        mol (list): Two-element list with the source of each molecule; each
            entry should match the corresponding ``mol_operation_type``.
        index (list): Two-element list of atom indices, one per molecule, at
            which the two molecules are expected to bind (i.e. the atoms
            brought into contact when linking them into the complex).
        bond_order (int, optional): Bond order to assume between the two
            linked atoms when building the complex; defaults to 1.
        db (str or dict, optional): Database credentials; path to db.json or a
            dict; if ``None``, read from the configuration files.
        name (str, optional): Name of the workflow; defaults to
            "binding_energy_calculation".
        working_dir (str, optional): Working directory for input/output files;
            defaults to the current working directory.
        opt_gaussian_inputs (dict, optional): Parameters for the optimization
            steps; defaults to B3LYP/6-31G(d).
        freq_gaussian_inputs (dict, optional): Parameters for the frequency
            steps; defaults to B3LYP/6-31G(d).
        solvent_gaussian_inputs (str, optional): Gaussian-style implicit
            solvent string (e.g. "(Solvent=Water)"); translated internally into
            ORCA's CPCM option, and also recorded as-is in the final db
            document for consistency with the Gaussian workflow's metadata.
        solvent_properties (dict, optional): Additional solvent properties
            recorded in the final db document (e.g. {"EPS": 12}).
        cart_coords (bool, optional): Accepted for interface parity; the ORCA
            backend only supports cartesian-coordinate inputs (``True``).
        oxidation_states (dict, optional): Oxidation states used to derive
            molecule charges (e.g. {"Li": 1, "O": -2}).
        skips (list, optional): Two-element list of jobs to skip per molecule
            (e.g. ``[["opt"], None]``); defaults to ``[None, None]`` (skip
            nothing).
        orca_cmd (str, optional): Path to the ORCA executable; falls back to
            the ORCA_CMD environment variable, then to "orca" on PATH. Must be
            an absolute path when ``num_cores`` > 1.
        num_cores (int, optional): Number of parallel ORCA processes per
            calculation; defaults to 1.
        memory (int, optional): Memory per core in MB; defaults to 4000.
        kwargs (keyword arguments): Additional kwargs passed to
            ``BindingEnergytoDB``/``Workflow`` (e.g. ``tag``).

    Returns:
        Workflow

    """
    working_dir = working_dir or os.getcwd()
    tag = kwargs.get("tag", "unknown")
    solvent = _to_orca_solvent(solvent_gaussian_inputs)

    opt_gaussian_inputs = opt_gaussian_inputs or {
        "functional": "b3lyp",
        "basis_set": "6-31G(d)",
        "route_parameters": {"Opt": None},
    }
    freq_gaussian_inputs = freq_gaussian_inputs or {
        "functional": "b3lyp",
        "basis_set": "6-31G(d)",
        "route_parameters": {"Freq": None},
    }

    skips = skips or [None, None]

    # engine-level settings shared by every ORCA step in this workflow
    orca_settings = {}
    if orca_cmd:
        orca_settings["orca_cmd"] = orca_cmd
    if num_cores:
        orca_settings["num_cores"] = num_cores
    if memory:
        orca_settings["memory"] = memory

    # Fireworks 1-4: optimize + frequency for each of the two molecules
    processed_mol_1 = process_mol(mol_operation_type[0], mol[0], db=db)
    processed_mol_2 = process_mol(mol_operation_type[1], mol[1], db=db)

    _, label_1, fws_1 = common_fw(
        mol=processed_mol_1,
        working_dir=working_dir,
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        mol_name="mol_1",
        gout_key="mol_1",
        skips=skips[0],
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        solvent=solvent,
        tag=tag,
        **orca_settings,
    )

    _, label_2, fws_2 = common_fw(
        mol=processed_mol_2,
        working_dir=working_dir,
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        mol_name="mol_2",
        gout_key="mol_2",
        skips=skips[1],
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        solvent=solvent,
        tag=tag,
        **orca_settings,
    )

    # Firework 5: combine the two optimized molecules and optimize the complex
    # (linking and optimizing must be one Firework -- see module docstring)
    linked_opt_fw = LinkedMolOrcaFW(
        gout_keys=["mol_1", "mol_2"],
        index=index,
        bond_order=bond_order,
        db=db,
        name="linked_mol_optimization",
        parents=fws_1 + fws_2,
        working_dir=os.path.join(working_dir, "linked_mol", "Optimization"),
        gaussian_input_params=opt_gaussian_inputs,
        gout_key="mol_linked_opt",
        cart_coords=cart_coords,
        solvent=solvent,
        tag=tag,
        **orca_settings,
    )

    # Firework 6: frequency calculation on the optimized complex
    linked_freq_fw = OrcaFW(
        prev_calc_key="mol_linked_opt",
        gaussian_input_params=freq_gaussian_inputs,
        db=db,
        name="linked_mol_frequency",
        parents=linked_opt_fw,
        working_dir=os.path.join(working_dir, "linked_mol", "Frequency"),
        gout_key="mol_linked",
        cart_coords=cart_coords,
        solvent=solvent,
        tag=tag,
        **orca_settings,
    )

    fws = fws_1 + fws_2 + [linked_opt_fw, linked_freq_fw]

    # Firework 7: gather energies, compute binding energy, save to db
    fw_analysis = Firework(
        BindingEnergytoDB(
            index=index,
            db=db,
            keys=["mol_1", "mol_2", "mol_linked"],
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            **{
                i: j
                for i, j in kwargs.items()
                if i
                in BindingEnergytoDB.required_params
                + BindingEnergytoDB.optional_params
            },
        ),
        parents=fws[:],
        name="binding_energy_analysis",
        spec={
            "tag": tag,
            "_launch_dir": os.path.join(working_dir, "analysis"),
        },
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name="{}_{}".format(name, "_".join([label_1, label_2])),
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS},
    )
