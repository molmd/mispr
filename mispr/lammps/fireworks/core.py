# coding: utf-8


# Defines core fireworks for ambertools and lammps.

import os
import logging

from fireworks import Firework

from mdproptools.structural.rdf_cn import calc_atomic_rdf, calc_molecular_rdf

from mispr.gaussian.workflows.base.core import _process_mol_check
from mispr.gaussian.firetasks.geo_transformation import ProcessMoleculeInput
from mispr.lammps.firetasks.run import RunTleap, RunLammpsDirect, RunParmchk, RunAntechamber
from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.lammps.firetasks.write_inputs import (
    LabelFFDict,
    WriteControlFile,
    WriteTleapScript,
    LabelFFDictFromDB,
)
from mispr.lammps.firetasks.parse_outputs import GetRDF, CalcDiff, ProcessPrmtop

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames


def ambertools_tasks(**kwargs):
    common_t = [
        RunAntechamber(
            **{
                i: j
                for i, j in kwargs.items()
                if i in RunAntechamber.required_params + RunAntechamber.optional_params
            }
        ),
        RunParmchk(
            **{
                i: j
                for i, j in kwargs.items()
                if i in RunParmchk.required_params + RunParmchk.optional_params
            }
        ),
        WriteTleapScript(
            **{
                i: j
                for i, j in kwargs.items()
                if i
                in WriteTleapScript.required_params + WriteTleapScript.optional_params
            }
        ),
        RunTleap(
            **{
                i: j
                for i, j in kwargs.items()
                if i in RunTleap.required_params + RunTleap.optional_params
            }
        ),
        ProcessPrmtop(
            **{
                i: j
                for i, j in kwargs.items()
                if i in ProcessPrmtop.required_params + ProcessPrmtop.optional_params
            }
        ),
    ]
    return common_t


class GetFFDictFW(Firework):
    def __init__(
        self,
        mol,
        mol_operation_type,
        data,
        operation_type="get_from_esp",
        label="",
        name="get_ff_dict",
        parents=None,
        working_dir=None,
        db=None,
        save_ff_to_db=False,
        save_ff_to_file=True,
        ff_filename="ff.json",
        tag="unknown",
        **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        # TODO: add logic for ensuring that directory name is legal or at \
        #  least usable
        os.makedirs(working_dir, exist_ok=True)

        tasks.append(
            ProcessMoleculeInput(
                mol=mol,
                operation_type=mol_operation_type,
                db=db,
                working_dir=working_dir,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in ProcessMoleculeInput.required_params
                    + ProcessMoleculeInput.optional_params
                }
            )
        )

        if operation_type == "get_from_esp":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the esp "
                    'file for operation_type="get_from_esp"'
                )
            tasks += ambertools_tasks(
                db=db,
                working_dir=working_dir,
                input_filename_a=data,
                unique_molecule_name=label,
                save_ff_to_db=save_ff_to_db,
                save_ff_to_file=save_ff_to_file,
                ff_filename=ff_filename,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in RunAntechamber.required_params
                    + RunAntechamber.optional_params
                    + RunParmchk.required_params
                    + RunParmchk.optional_params
                    + WriteTleapScript.required_params
                    + WriteTleapScript.optional_params
                    + RunTleap.required_params
                    + RunTleap.optional_params
                    + ProcessPrmtop.required_params
                    + ProcessPrmtop.optional_params
                }
            )

        elif operation_type == "get_from_prmtop":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the prmtop "
                    'file for operation_type="get_from_prmtop"'
                )
            tasks.append(
                ProcessPrmtop(
                    working_dir=working_dir,
                    db=db,
                    prmtop_path=data,
                    unique_molecule_name=label,
                    save_ff_to_db=save_ff_to_db,
                    save_ff_to_file=save_ff_to_file,
                    ff_filename=ff_filename,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in ProcessPrmtop.required_params + ProcessPrmtop.optional_params
                    }
                )
            )

        elif operation_type == "get_from_dict":
            if not isinstance(data, dict):
                raise TypeError(
                    "data must be a dict for " 'operation_type="get_from_dict"'
                )
            tasks.append(
                LabelFFDict(
                    unlabeled_dict=data,
                    working_dir=working_dir,
                    label=label,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in LabelFFDict.required_params + LabelFFDict.optional_params
                    }
                )
            )

        elif operation_type == "get_from_file":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the ff "
                    'json file for operation_type="get_from_file"'
                )
            tasks.append(
                LabelFFDict(
                    ff_file=data,
                    working_dir=working_dir,
                    label=label,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in LabelFFDict.required_params + LabelFFDict.optional_params
                    }
                )
            )

        elif operation_type == "get_from_db":
            tasks.append(
                LabelFFDictFromDB(
                    filter=data,
                    db=db,
                    working_dir=working_dir,
                    label=label,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in LabelFFDictFromDB.required_params
                        + LabelFFDictFromDB.optional_params
                    }
                )
            )

        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(GetFFDictFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )


class RunLammpsFW(Firework):
    def __init__(
        self,
        control_file=None,
        db=None,
        name="run_lammps",
        parents=None,
        working_dir=None,
        save_run_to_db=True,
        save_run_to_file=False,
        tag="unknown",
        **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        if not control_file:
            tasks.append(
                WriteControlFile(
                    working_dir=working_dir,
                    db=db,
                    save_runs_to_db=save_run_to_db,
                    save_runs_to_file=save_run_to_file,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in WriteControlFile.required_params
                        + WriteControlFile.optional_params
                    }
                )
            )

        tasks.append(
            RunLammpsDirect(
                working_dir=working_dir,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i in RunLammpsDirect.required_params + RunLammpsDirect.optional_params
                }
            )
        )

        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})

        super(RunLammpsFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )


class RunAnalysisFW(Firework):
    def __init__(
        self,
        md_property,
        name="run_analysis",
        parents=None,
        working_dir=None,
        tag="unknown",
        **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        if md_property == "diffusion":
            diff_kwargs = {
                i: j
                for i, j in kwargs.items()
                if i in CalcDiff.required_params + CalcDiff.optional_params
            }
            diff_kwargs.update({"working_dir": working_dir})
            tasks.append(CalcDiff(**diff_kwargs))

        elif md_property == "rdf":
            rdf_kwargs = {
                i: j
                for i, j in kwargs.items()
                if i
                in inspect.getfullargspec(calc_atomic_rdf).args
                + inspect.getfullargspec(calc_molecular_rdf).args
                + GetRDF.required_params
                + GetRDF.optional_params
            }
            rdf_kwargs.update({"working_dir": working_dir})
            tasks.append(GetRDF(**rdf_kwargs))
        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(RunAnalysisFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )
