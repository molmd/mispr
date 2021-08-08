# coding: utf-8


# Defines custom fireworks for ambertools and lammps.

import os
import inspect
import logging

from fireworks import Firework

from mdproptools.structural.rdf_cn import calc_atomic_rdf, calc_molecular_rdf

from mispr.lammps.firetasks.run import RunTleap, RunLammps, RunParmchk, RunAntechamber
from mispr.gaussian.utilities.metadata import get_mol_formula
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


class EspToData(Firework):
    def __init__(
        self, prmtop_filename, file_label, working_dir=None,
    ):
        pass


class GetFFDictFW(Firework):
    def __init__(
        self,
        mol,
        data,
        operation_type="get_from_esp",
        db=None,
        label="",
        name="get_ff_dict",
        parents=None,
        working_dir=None,
        **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        if not label:
            label = get_mol_formula(mol)

        if operation_type == "get_from_esp":
            if not isinstance(data, str):
                raise TypeError(
                    "data must be a str of the path to the esp file"
                    'for operation_type="get_from_esp"'
                )
            tasks += ambertools_tasks(
                mol,
                db=db,
                working_dir=working_dir,
                input_filename_a=data,
                unique_molecule_name=label,
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
                    "data must be a str of the path to the prmtop file"
                    'for operation_type="get_from_prmtop"'
                )
            tasks.append(
                ProcessPrmtop(
                    molecule=mol,
                    working_dir=working_dir,
                    db=db,
                    prmtop_path=data,
                    unique_molecule_name=label,
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
                    'data must be a dict for \
                                operation_type="get_from_dict"'
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


class WriteLammpsInputs(Firework):
    def __init__(
        self,
        system_mixture_data_type,
        system_mixture_data,
        system_box_data,
        system_box_data_type,
        control_template,
        control_settings,
        working_dir=None,
        data_filename="complex.data",
        control_file="control.lammpsin",
        **kwargs
    ):
        """

        :param system_mol_data: [dict]
                      {
                      mol1_label: {"molecule": pmg.Molecule,
                                   "ff_param_method": str (default: "get_from_esp"),
                                   "ff_param_data": str or dict,
                                   "mol_mixture_type": "Solutes" or "Solvents",
                                   "mixture_data": int or dict},
                      ...,
                      moln_label: {...}
                      }

        :param ff_gen_type:
        :param working_dir:
        :param data_file:
        :param control_file:
        """
        tasks = []
        working_dir = working_dir or os.getcwd()

        pass


class RunLammpsFW(Firework):
    def __init__(
        self,
        control_file=None,
        db=None,
        name="run_lammps",
        parents=None,
        working_dir=None,
        **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        if not control_file:
            tasks.append(
                WriteControlFile(
                    working_dir=working_dir,
                    db=db,
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
        self, property, name="run_analysis", parents=None, working_dir=None, **kwargs
    ):
        tasks = []
        working_dir = working_dir or os.getcwd()

        if property == "diffusion":
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

        elif property == "rdf":
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


if __name__ == "__main__":

    from pymatgen.io.gaussian import GaussianOutput
    from fireworks import FWorker, LaunchPad
    from fireworks.core.rocket_launcher import launch_rocket

    # wd = os.getcwd()
    # os_path = os.path.normpath(os.path.join(wd, 'test', 'test2'))
    # print(os_path)
    # pathlib_path = pathlib.Path(os_path)
    # print(pathlib_path)
    # pathlib_path.mkdir(parents = True, exist_ok = True)

    working_dir = "/Users/matt/Documents/GitHub/infrastructure/infrastructure/lammps/tests/test_files/antechamber"
    esp_file_path = "/Users/matt/Documents/GitHub/infrastructure/infrastructure/lammps/tests/test_files/antechamber/dhps.esp"
    gout_file_path = "/Users/matt/Documents/GitHub/infrastructure/infrastructure/lammps/tests/test_files/dhps.out"
    mol2_file = "mol.mol2"

    dhps_gout = GaussianOutput(gout_file_path)
    dhps_mol = dhps_gout.structures[-1]
    dhps_mol.set_charge_and_spin(dhps_gout.charge, dhps_gout.spin_multiplicity)

    launchpad = LaunchPad(
        host="mongodb+srv://mbliss01:idlewide@gettingstarted.dt0sv.mongodb.net/fireworks",
        uri_mode=True,
    )
    launchpad.reset("", require_password=False)

    # firework = GetFFDictFW(dhps_mol,
    #                        esp_file_path,
    #                        working_dir=working_dir,
    #                        output_filename_a=mol2_file,
    #                        prmtop_filename=mol2_file.split('.')[0] + ".prmtop")

    working_dir = "/Users/matt/Documents/GitHub/infrastructure/infrastructure/lammps/tests/test_files/analysis"
    log_file_name = "log.lammps"
    firework = RunAnalysisFW(
        "diffusion",
        msd_method="from_log",
        file_pattern=log_file_name,
        working_dir=working_dir,
        kwargs={"dt": 2},
    )

    # dump_filename = "dump.0.1M_PHEN_1.3_NaOH_SPCE.npt.dump"
    # r_cut = 10
    # bin_size = 0.1
    # num_types = 8
    # mass = [14.01, 12.01, 1.008, 16.0, 1.008, 16.0, 1.008, 22.99]
    # partial_relations = [[1, 1], [1, 8]]
    # csv_filename = "rdf.csv"
    # firework = RunAnalysisFW("rdf",
    #                          working_dir = working_dir,
    #                          rdf_type = "atomic",
    #                          kwargs = {"filename": os.path.join(working_dir, dump_filename),
    #                                    "r_cut": r_cut,
    #                                    "bin_size": bin_size,
    #                                    "num_types": num_types,
    #                                    "mass": mass,
    #                                    "partial_relations": partial_relations})

    launchpad.add_wf(firework)
    launch_rocket(launchpad, FWorker())
