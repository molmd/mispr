# coding: utf-8


# Defines common fireworks used in Gaussian workflows.

import os
import logging

from fireworks import Firework

from mispr.gaussian.firetasks.run_calc import RunGaussianCustodian
from mispr.gaussian.firetasks.write_inputs import WriteInput
from mispr.gaussian.firetasks.parse_outputs import ProcessRun, RetrieveGaussianOutput
from mispr.gaussian.firetasks.geo_transformation import ProcessMoleculeInput

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.3"

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames


def common_tasks(
    db,
    input_file,
    output_file,
    gaussian_input_params,
    cart_coords,
    oxidation_states,
    **kwargs
):
    """
    Defines a list of common tasks for Gaussian fireworks, i.e.
    writing an input file, running the calculation, and
    parsing the output.

    Args:
        db (str or dict): database credentials to store the run; could
            be provided as the path to the db.json file or in the form
            of a dictionary
        input_file (str): name of the input file to be written
        output_file (str): name of the Gaussian output file
        gaussian_input_params (dict): a dictionary of parameters to be
            used in creating the Gaussian input file
        cart_coords (bool): whether to write cartesian coordinates or
            not; default is True
        oxidation_states (dict): a dictionary of element symbols and
            their oxidation states used in setting the charge on the
            molecule
        **kwargs: other kwargs that are passed to:
            1. mispr.gaussian.firetasks.write_inputs.WriteInput
            2. mispr.gaussian.firetasks.run_calc.RunGaussianCustodian
            3. mispr.gaussian.firetasks.parse_outputs.ProcessRun

    Returns:
        list of Firetasks
    """
    return [
        WriteInput(
            input_file=input_file,
            gaussian_input_params=gaussian_input_params,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            **{
                i: j
                for i, j in kwargs.items()
                if i in WriteInput.required_params + WriteInput.optional_params
            }
        ),
        RunGaussianCustodian(
            input_file=input_file,
            output_file=output_file,
            **{
                i: j
                for i, j in kwargs.items()
                if i
                in RunGaussianCustodian.required_params
                + RunGaussianCustodian.optional_params
            }
        ),
        ProcessRun(
            run=output_file,
            operation_type="get_from_gout_file",
            db=db,
            input_file=input_file,
            **{
                i: j
                for i, j in kwargs.items()
                if i in ProcessRun.required_params + ProcessRun.optional_params
            }
        ),
    ]


class CalcFromMolFW(Firework):
    """
    Runs a Gaussian calculation from a molecule.
    """
    def __init__(
        self,
        mol,
        mol_operation_type="get_from_mol",
        db=None,
        name="calc_from_mol",
        parents=None,
        working_dir=None,
        input_file="mol.com",
        output_file="mol.out",
        gaussian_input_params={},
        cart_coords=True,
        oxidation_states=None,
        tag="unknown",
        **kwargs
    ):
        """
        Args:
            mol (Molecule, GaussianOutput, str, dict): source of the
                molecule to be processed. Should match the mol_operation_type
            mol_operation_type (str): the type of molecule operation.
                See process_mol defined in mispr/gaussian/utilities/mol.py
                for supported operations; defaults to "get_from_mol"
            db (str or dict): database credentials; could be provided as
                the path to the db.json file or in the form of a dictionary;
                if none is provided, attempts to get it from the
                configuration files
            name (str): name of the Firework; defaults to "calc_from_mol"
            parents (Firework or [Firework]): list of parent FWs this FW
                depends on
            working_dir (str): working directory for the calculation;
                defaults to the current directory
            input_file (str): name of the Gaussian input file to be created;
                defaults to "mol.com"
            output_file (str): name of the Gaussian output file to be output;
                defaults to "mol.out"
            gaussian_input_params (dict): dictionary of parameters to be
                used in the Gaussian input file
            cart_coords (bool): whether the coordinates are cartesian or
                z-matrix; defaults to True
            oxidation_states (list): list of oxidation states for each
                atom; defaults to None
            tag (str): tag for the calculation; the provided tag will be
                stored in the db documents for easy retrieval; defaults
                to "unknown"
            **kwargs: other kwargs that are passed to:
                1. Firework.__init__.
                2. mispr.gaussian.firetasks.geo_transformation.ProcessMoleculeInput
                3. mispr.gaussian.fireworks.common_tasks
        """

        t = []
        working_dir = working_dir or os.getcwd()
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        t.append(
            ProcessMoleculeInput(
                mol=mol,
                operation_type=mol_operation_type,
                db=db,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in ProcessMoleculeInput.required_params
                    + ProcessMoleculeInput.optional_params
                }
            )
        )
        t += common_tasks(
            db,
            input_file,
            output_file,
            gaussian_input_params,
            cart_coords,
            oxidation_states,
            # additional_fw=CalcFromMolFW(mol=mol,
            #                             mol_operation_type=mol_operation_type,
            #                             db=db,
            #                             name=name,
            #                             parents=parents,
            #                             working_dir=working_dir,
            #                             input_file=input_file,
            #                             output_file=output_file,
            #                             gaussian_input_params=gaussian_input_params,
            #                             cart_coords=cart_coords,
            #                             oxidation_states=oxidation_states,
            #                             tag=tag,
            #                             **kwargs),
            **kwargs
        )
        spec = kwargs.pop("spec", {})
        spec.update(
            {"tag": tag, "_launch_dir": working_dir, "_add_launchpad_and_fw_id": True}
        )
        super(CalcFromMolFW, self).__init__(
            t,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )


class CalcFromRunsDBFW(Firework):
    """
    Runs a Gaussian calculation from a previous calculation or the
    runs database.
    """
    def __init__(
        self,
        db=None,
        name="calc_from_runs_db",
        parents=None,
        gaussian_input_params=None,
        working_dir=None,
        input_file="mol.com",
        output_file="mol.out",
        cart_coords=True,
        tag="unknown",
        **kwargs
    ):
        """
        Args:
            db (str or dict): database credentials; could be provided as
                the path to the db.json file or in the form of a dictionary;
                if none is provided, attempts to get it from the
                configuration files
            name (str): name of the Firework; defaults to "calc_from_runs_db"
            parents (Firework or [Firework]): list of parent FWs this FW
                depends on
            gaussian_input_params (dict): dictionary of parameters to be
                used in the Gaussian input file
            working_dir (str): working directory for the calculation;
                defaults to the current directory
            input_file (str): name of the Gaussian input file to be created;
                defaults to "mol.com"
            output_file (str): name of the Gaussian output file to be output;
                defaults to "mol.out"
            cart_coords (bool): whether the coordinates are cartesian or
                z-matrix; defaults to True
            tag (str): tag for the calculation; the provided tag will be
                stored in the db documents for easy retrieval; defaults
                to "unknown"
            **kwargs: other kwargs that are passed to:
                1. Firework.__init__.
                2. mispr.gaussian.firetasks.parse_outputs.RetrieveGaussianOutput
                3. mispr.gaussian.fireworks.common_tasks
        """
        t = []
        working_dir = working_dir or os.getcwd()
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        t.append(
            RetrieveGaussianOutput(
                db=db,
                gaussian_input_params=gaussian_input_params,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in RetrieveGaussianOutput.required_params
                    + RetrieveGaussianOutput.optional_params
                }
            )
        )
        t += common_tasks(
            db,
            input_file,
            output_file,
            gaussian_input_params,
            cart_coords,
            oxidation_states=None,
            **kwargs
        )
        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(CalcFromRunsDBFW, self).__init__(
            t,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )
