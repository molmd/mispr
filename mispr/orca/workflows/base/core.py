"""Define a list of common Fireworks used in ORCA workflows."""

import os
import logging

from fireworks import Workflow

from mispr.orca.fireworks.core import OrcaFW

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"

logger = logging.getLogger(__name__)

WORKFLOW_KWARGS = Workflow.__init__.__code__.co_varnames


def common_fw(
    mol=None,
    prev_calc_key=None,
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    gout_key=None,
    db=None,
    mol_name=None,
    dir_head=None,
    dir_structure=None,
    skips=None,
    tag="unknown",
    cart_coords=True,
    oxidation_states=None,
    solvent=None,
    **kwargs,
):
    """Define a list of Fireworks commonly used in ORCA workflows: optimize a molecule, then run a frequency analysis on the optimized structure.

    Mirrors the Gaussian original, swapping the engine Firework for
    ``OrcaFW``.

    Args:
        mol (Molecule, optional): pymatgen Molecule to run the calculations on;
            required unless ``prev_calc_key`` is given (e.g. when chaining onto
            a state produced by an earlier node in a tree-structured workflow,
            such as ``mispr.orca.workflows.base.ip_ea``, where the geometry
            only exists once that earlier Firework actually runs).
        prev_calc_key (str, optional): Key of a previous run in
            fw_spec["gaussian_output"] to use as the starting structure instead
            of ``mol``.
        working_dir (str): Path of the working directory.
        opt_gaussian_inputs (dict): Parameters for the optimization step (see
            ``mispr.orca.fireworks.core.OrcaFW``).
        freq_gaussian_inputs (dict): Parameters for the frequency step.
        gout_key (str, optional): Unique key for the run dict; defaults to
            "mol".
        db (str or dict, optional): Database credentials.
        mol_name (str, optional): Name of the molecule, used for labeling.
        dir_head (str, optional): Head of the directory where the workflow will
            run; defaults to ``mol_name``.
        dir_structure (list, optional): Additional subfolders appended under
            ``dir_head``.
        skips (list, optional): List of jobs to skip; e.g. ["opt"] or ["freq"].
        tag (str, optional): Tag stored in the db documents.
        cart_coords (bool, optional): Accepted for interface parity; the ORCA
            backend only supports cartesian-coordinate inputs (``True``).
        oxidation_states (dict, optional): Oxidation states used to derive the
            molecule's charge in the optimization step (e.g. {"Li": 1, "O": -2}).
        solvent (dict, optional): Implicit solvent options forwarded to
            ``OrcaFW`` (e.g. {"solvent": "water"}); if not provided, gas phase.
        kwargs: Additional kwargs passed through to ``OrcaFW`` (e.g.
            ``orca_cmd``, ``num_cores``); unrecognized ones are ignored
            (accepted for interface parity with the Gaussian ``common_fw``).

    Returns:
        Molecule, str, list: The input molecule, a label, and a list of
        Fireworks.

    """
    fws = []
    if not gout_key:
        gout_key = "mol"

    label = mol_name or "mol"
    dir_head = dir_head or label
    dir_structure = dir_structure or []
    working_dir = os.path.join(working_dir, dir_head, *dir_structure)

    skips = skips or []

    if not skips:
        opt_fw = OrcaFW(
            molecule=mol,
            prev_calc_key=prev_calc_key,
            gaussian_input_params=opt_gaussian_inputs,
            db=db,
            name=f"{label}_optimization",
            working_dir=os.path.join(working_dir, "Optimization"),
            gout_key=gout_key + "_opt",
            tag=tag,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            solvent=solvent,
            **kwargs,
        )
        fws.append(opt_fw)

        freq_fw = OrcaFW(
            prev_calc_key=gout_key + "_opt",
            gaussian_input_params=freq_gaussian_inputs,
            db=db,
            name=f"{label}_frequency",
            parents=opt_fw,
            working_dir=os.path.join(working_dir, "Frequency"),
            gout_key=gout_key,
            tag=tag,
            cart_coords=cart_coords,
            solvent=solvent,
            **kwargs,
        )
        fws.append(freq_fw)

    elif len(skips) == 1 and skips[0].lower() == "opt":
        freq_fw = OrcaFW(
            molecule=mol,
            prev_calc_key=prev_calc_key,
            gaussian_input_params=freq_gaussian_inputs,
            db=db,
            name=f"{label}_frequency",
            working_dir=os.path.join(working_dir, "Frequency"),
            gout_key=gout_key,
            tag=tag,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            solvent=solvent,
            **kwargs,
        )
        fws.append(freq_fw)

    elif len(skips) == 1 and skips[0].lower() == "freq":
        opt_fw = OrcaFW(
            molecule=mol,
            prev_calc_key=prev_calc_key,
            gaussian_input_params=opt_gaussian_inputs,
            db=db,
            name=f"{label}_optimization",
            working_dir=os.path.join(working_dir, "Optimization"),
            gout_key=gout_key + "_opt",
            tag=tag,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            solvent=solvent,
            **kwargs,
        )
        fws.append(opt_fw)

    return mol, label, fws
