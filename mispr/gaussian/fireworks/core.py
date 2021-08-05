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
__version__ = "0.0.1"

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
