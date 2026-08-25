"""Define the electrostatic partial charge (ESP/CHELPG) workflow, using ORCA instead of Gaussian.

Mirrors ``mispr.gaussian.workflows.base.esp.get_esp_charges``: optimize +
frequency the molecule, run a CHELPG ESP single-point on the optimized
geometry, then save the result.
"""

import os

from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.firetasks.parse_outputs import ESPtoDB
from mispr.orca.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.orca.firetasks.run_calc import ESP
from fireworks import Firework, Workflow

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"


def get_esp_charges(
    mol_operation_type,
    mol,
    db=None,
    name="esp_charges_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    esp_gaussian_inputs=None,
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
    """Define a workflow for calculating the electrostatic partial charges with ORCA.

    * **Firework 1**: Optimize the molecule.
    * **Firework 2**: Run a frequency analysis.
    * **Firework 3**: Run a CHELPG ESP calculation.
    * **Firework 4**: Create ESP document/json file.

    Args:
        mol_operation_type (str): The type of molecule operation. See
            ``process_mol`` defined in ``mispr/gaussian/utilities/mol.py`` for
            supported operations.
        mol (Molecule, GaussianOutput, str, dict): Source of the molecule to be
            processed. Should match the ``mol_operation_type``.
        db (str or dict, optional): Database credentials; could be provided as
            the path to the "db.json" file or in the form of a dictionary; if
            none is provided, attempts to get it from the configuration files.
        name (str, optional): Name of the workflow. Defaults to
            "esp_charges_calculation".
        working_dir (str, optional): Path of the working directory where any
            required input files can be found and output will be created.
            Defaults to the current working directory.
        opt_gaussian_inputs (dict, optional): Input parameters for the
            optimization step, in the shared Gaussian-style format (see
            ``mispr.orca.fireworks.core.OrcaFW``); basis set names must be ones
            ORCA recognizes. Defaults to B3LYP/6-31G(d).
        freq_gaussian_inputs (dict, optional): Input parameters for the
            frequency step; default parameters will be used if not specified.
        esp_gaussian_inputs (dict, optional): Input parameters for the ESP step;
            recognized keys are "functional" (used as ``method_esp``) and
            "basis_set" (used as ``basis_esp``). Defaults to HF/6-31G*.
        solvent_gaussian_inputs (str, optional): Accepted for interface parity
            with the Gaussian workflow; currently unused by the ORCA backend
            (pass a ``solvent`` dict through ``kwargs`` instead to request
            CPCM).
        solvent_properties (dict, optional): Accepted for interface parity;
            currently unused by the ORCA backend.
        cart_coords (bool, optional): Accepted for interface parity; the ORCA
            backend only supports cartesian-coordinate inputs (``True``).
        oxidation_states (dict, optional): Dictionary of oxidation states used
            in setting the charge and spin multiplicity of the molecule; e.g.:
            {"Li": 1, "O": -2}. Defaults to None.
        skips (list, optional): List of jobs to skip; e.g.: ["opt", "freq"];
            defaults to None.
        orca_cmd (str, optional): Path to the ORCA executable; falls back to
            the ORCA_CMD environment variable, then to "orca" on PATH. Must be
            an absolute path when ``num_cores`` > 1.
        num_cores (int, optional): Number of parallel ORCA processes; defaults
            to 1.
        memory (int, optional): Memory per core in MB; defaults to 4000.
        kwargs (keyword arguments): Additional kwargs to be passed to the
            workflow.

    Returns:
        tuple:
            - Workflow
            - str: Label of the molecule (e.g. "H2O", "water", etc.).

    """
    working_dir = working_dir or os.getcwd()
    processed_mol = process_mol(mol_operation_type, mol, db=db)

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
    esp_gaussian_inputs = esp_gaussian_inputs or {}

    tag = kwargs.get("tag", "unknown")

    # engine-level settings shared by every ORCA step in this workflow
    orca_settings = {}
    if orca_cmd:
        orca_settings["orca_cmd"] = orca_cmd
    if num_cores:
        orca_settings["num_cores"] = num_cores
    if memory:
        orca_settings["memory"] = memory

    _, label, fws = common_fw(
        mol=processed_mol,
        working_dir=working_dir,
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        mol_name="mol",
        gout_key="mol",
        skips=skips,
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        solvent=None,
        tag=tag,
        **orca_settings,
    )

    # spec must carry "tag" itself -- unlike the opt/freq Fireworks built by
    # common_fw (which get it via OrcaFW's spec.update({"tag": tag, ...})), a
    # plain Firework doesn't inherit spec from its parents, so ESPtoDB would
    # otherwise KeyError on fw_spec["tag"] downstream
    esp_task_kwargs = dict(orca_settings)
    if "functional" in esp_gaussian_inputs:
        esp_task_kwargs["method_esp"] = esp_gaussian_inputs["functional"]
    if "basis_set" in esp_gaussian_inputs:
        esp_task_kwargs["basis_esp"] = esp_gaussian_inputs["basis_set"]
    esp_working_dir = os.path.join(working_dir, "mol", "ESP")
    if not os.path.exists(esp_working_dir):
        os.makedirs(esp_working_dir)
    esp_fw = Firework(
        tasks=[
            ESP(prev_calc_key="mol", gout_key="mol_esp", db=db, **esp_task_kwargs)
        ],
        parents=fws[:],
        name=f"{label}_esp",
        spec={"tag": tag, "_launch_dir": esp_working_dir},
    )
    fws.append(esp_fw)

    fw_analysis = Firework(
        ESPtoDB(
            db=db,
            keys=["mol", "mol_esp"],
            **{
                i: j
                for i, j in kwargs.items()
                if i in ESPtoDB.required_params + ESPtoDB.optional_params
            },
        ),
        parents=fws[:],
        name=f"{label}_esp_analysis",
        spec={"tag": tag},
    )
    fws.append(fw_analysis)

    return (
        Workflow(
            fws,
            name=name,
            **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS},
        ),
        label,
    )
