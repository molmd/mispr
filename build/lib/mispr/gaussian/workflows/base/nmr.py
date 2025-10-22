"""Define the nuclear magnetic resonance workflow."""

import os
import logging

from fireworks import Firework, Workflow

from mispr.gaussian.fireworks.core import CalcFromRunsDBFW
from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.utilities.inputs import handle_gaussian_inputs
from mispr.gaussian.utilities.metadata import get_job_name
from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.gaussian.firetasks.parse_outputs import NMRtoDB

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.4"


logger = logging.getLogger(__name__)


def get_nmr_tensors(
    mol_operation_type,
    mol,
    db=None,
    name="nmr_tensor_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    nmr_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    **kwargs
):
    """
    Define a workflow for calculating the nuclear magnetic resonance tensors.

    * **Firework 1**: Optimize the molecule.
    * **Firework 2**: Calculate the frequencies.
    * **Firework 3**: Calculate the NMR tensors.
    * **Firework 4**: Create NMR document/json file.

    Args:
        mol_operation_type (str): The type of molecule operation. See ``process_mol``
            defined in ``mispr/gaussian/utilities/mol.py`` for supported operations.
        mol (Molecule, GaussianOutput, str, dict): Source of the molecule to be
            processed. Should match the ``mol_operation_type``.
        db (str or dict, optional): Database credentials; could be provided as the path
            to the "db.json" file or in the form of a dictionary; if none is provided,
            attempts to get it from the configuration files.
        name (str, optional): Name of the workflow. Defaults to "nmr_tensor_calculation".
        working_dir (str, optional): Path of the working directory where any required
            input files can be found and output will be created. Defaults to the current
            working directory.
        opt_gaussian_inputs (dict, optional): Dictionary of Gaussian input parameters;
            e.g.:

            .. code-block:: python

                {
                    "functional": "B3LYP",
                    "basis_set": "6-31G(d)",
                    "route_parameters": {"Opt": None},
                    "link0_parameters": {
                        "%chk": "checkpoint.chk",
                        "%mem": "45GB",
                        "%NProcShared": "24"}
                }

            The above default parameters will be used if not specified.
        freq_gaussian_inputs (dict, optional): Dictionary of Gaussian input parameters
            for the frequency step; default parameters will be used if not specified.
        nmr_gaussian_inputs (dict, optional): Dictionary of Gaussian input parameters
            for the NMR step; default parameters will be used if not specified.
        solvent_gaussian_inputs (str, optional): Gaussian input parameters corresponding
            to the implicit solvent model to be used in the NMR calculations, if any;
            e.g.:

            .. code-block:: python

                "(Solvent=TetraHydroFuran)"

            These parameters should only be specified here and not included in the main
            gaussian_inputs dictionary for each job (i.e. ``opt_gaussian_inputs``,
            ``freq_gaussian_inputs``, etc.). Defaults to None.
        solvent_properties (dict, optional): Additional input parameters to be used in
            the NMR calculations and relevant to the solvent model, if any; e.g.,
            {"EPS":12}. Defaults to None.
        cart_coords (bool, optional): Uses cartesian coordinates in writing Gaussian
            input files if set to ``True``, otherwise uses z-matrix. Defaults to ``True``.
        oxidation_states (dict, optional): Dictionary of oxidation states that can be
            used in setting the charge and spin multiplicity of the molecule; e.g.:
            {"Li":1, "O":-2}. Defaults to None.
        skips (list, optional): List of jobs to skip; e.g.: ["opt", "freq"]; defaults
            to None.
        kwargs (keyword arguments): Additional kwargs to be passed to the workflow.

    Returns:
        Workflow
    """
    fws = []
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)
    gout_keys = ["mol", "mol_nmr"]

    gaussian_inputs = handle_gaussian_inputs(
        {
            "opt": opt_gaussian_inputs,
            "freq": freq_gaussian_inputs,
            "nmr": nmr_gaussian_inputs,
        },
        solvent_gaussian_inputs,
        solvent_properties,
    )
    opt_gaussian_inputs = gaussian_inputs["opt"]
    freq_gaussian_inputs = gaussian_inputs["freq"]
    nmr_gaussian_inputs = gaussian_inputs["nmr"]

    _, label, opt_freq_fws = common_fw(
        mol_operation_type=mol_operation_type,
        mol=mol,
        working_dir=working_dir,
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        gout_key=gout_keys[0],
        skips=skips,
        **kwargs
    )
    fws += opt_freq_fws

    spec = kwargs.pop("spec", {})
    if not skips or len(skips) == 1 and skips[0].lower() == "opt":
        spec.update(
            {"proceed": {"has_gaussian_completed": True, "stationary_type": "Minimum"}}
        )
    else:
        spec.update({"proceed": {"has_gaussian_completed": True}})

    nmr_fw = CalcFromRunsDBFW(
        db,
        input_file="{}_nmr.com".format(label),
        output_file="{}_nmr.out".format(label),
        name=get_job_name(label, "nmr"),
        parents=fws[:],
        gaussian_input_params=nmr_gaussian_inputs,
        working_dir=os.path.join(working_dir, label, "NMR"),
        cart_coords=cart_coords,
        gout_key=gout_keys[1],
        spec=spec,
        **kwargs
    )
    fws.append(nmr_fw)

    fw_analysis = Firework(
        NMRtoDB(
            db=db,
            keys=gout_keys,
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            **{
                i: j
                for i, j in kwargs.items()
                if i in NMRtoDB.required_params + NMRtoDB.optional_params
            }
        ),
        parents=fws[:],
        name="{}-{}".format(label, "nmr_analysis"),
        spec={"_launch_dir": os.path.join(working_dir, label, "analysis")},
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name=get_job_name(label, name),
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS}
    )
