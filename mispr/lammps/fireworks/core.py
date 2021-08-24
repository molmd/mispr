# coding: utf-8


# Defines core fireworks for ambertools and lammps.

import os
import logging

import fireworks as Firework

from mdproptools.structural.rdf_cn import calc_atomic_rdf, calc_molecular_rdf

from mispr.lammps.firetasks.run import RunTleap, RunLammps, RunParmchk, RunAntechamber
from mispr.gaussian.utilities.metadata import get_chem_schema
from mispr.lammps.firetasks.write_inputs import (
    LabelFFDict,
    WriteControlFile,
    WriteTleapScript,
    LabelFFDictFromDB,
)
from mispr.lammps.firetasks.parse_outputs import GetMSD, GetRDF, CalcDiff, ProcessPrmtop

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames


def ambertools_tasks(mol, **kwargs):
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
            molecule=mol,
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
        **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        # TODO: add logic for ensuring that directory name is legal or at \
        #  least usable
        os.makedirs(working_dir, exist_ok=True)

        if not label:
            label = get_chem_schema(mol)["formula_alphabetical"]

        if operation_type == "get_from_esp":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the esp "
                    'file for operation_type="get_from_esp"'
                )
            tasks += ambertools_tasks(
                mol,
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
                    molecule=mol,
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
                    mol=mol,
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

        super(GetFFDictFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
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
        **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        if not control_file:
            tasks.append(
                WriteControlFile(
                    working_dir=working_dir,
                    db=db,
                    save_to_db=save_run_to_db,
                    save_to_file=save_run_to_file,
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
            RunLammps(
                working_dir=working_dir,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i in RunLammps.required_params + RunLammps.optional_params
                }
            )
        )

        super(RunLammpsFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )


class RunAnalysisFW(Firework):
    def __init__(
        self, md_property, name="run_analysis", parents=None, working_dir=None, **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        if md_property == "diffusion":
            msd_kwargs = {
                i: j
                for i, j in kwargs.items()
                if i in GetMSD.required_params + GetMSD.optional_params
            }
            msd_kwargs.update({"working_dir": working_dir})
            tasks.append(GetMSD(**msd_kwargs))
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

        super(RunAnalysisFW, self).__init__(
            tasks,
            parents=parents,
            name=name,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )
