# coding: utf-8

# Defines custom fireworks for ambertools and lammps

import os
import pathlib
import logging
import inspect

import fireworks as fw
import pymatgen.core.structure as pmgcs
import infrastructure.lammps.firetasks.write_inputs as ilftwi
import infrastructure.lammps.firetasks.run as ilftr
import infrastructure.lammps.firetasks.parse_outputs as ilfpo
import analysis.lammps.structural.rdf_cn as alsrc
import infrastructure.gaussian.utils.utils as iguu

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = fw.Firework.__init__.__code__.co_varnames


def ambertools_tasks(mol, **kwargs):

    common_t = [ilftr.RunAntechamber(
                    **{i: j for i, j in kwargs.items() if i in
                       ilftr.RunAntechamber.required_params
                       + ilftr.RunAntechamber.optional_params}),
                ilftr.RunParmchk(
                    **{i: j for i, j in kwargs.items() if i in
                       ilftr.RunParmchk.required_params
                       + ilftr.RunParmchk.optional_params}),
                ilftwi.WriteTleapScript(
                    **{i: j for i, j in kwargs.items() if i in
                       ilftwi.WriteTleapScript.required_params
                       + ilftwi.WriteTleapScript.optional_params}),
                ilftr.RunTleap(
                    **{i: j for i, j in kwargs.items() if i in
                       ilftr.RunTleap.required_params
                       + ilftr.RunTleap.optional_params}),
                ilfpo.ProcessPrmtop(
                    molecule=mol,
                    **{i: j for i, j in kwargs.items() if i in
                       ilfpo.ProcessPrmtop.required_params
                       + ilfpo.ProcessPrmtop.optional_params})]
    return common_t


class GetFFDictFW(fw.Firework):
    def __init__(self,
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
                 **kwargs):
        tasks = []
        working_dir = working_dir or os.getcwd()

        # TODO: add logic for ensuring that directory name is legal or at \
        #  least usable
        os.makedirs(working_dir, exist_ok=True)

        if not label:
            label = iguu.get_chem_schema(mol)["formula_alphabetical"]

        # print(operation_type)

        if operation_type == "get_from_esp":
            if not isinstance(data, str):
                raise TypeError('data must be a str of the path to the esp ' 
                                'file for operation_type="get_from_esp"')
            tasks += ambertools_tasks(
                          mol,
                          db=db,
                          working_dir=working_dir,
                          input_filename_a=data,
                          unique_molecule_name=label,
                          save_ff_to_db=save_ff_to_db,
                          save_ff_to_file=save_ff_to_file,
                          ff_filename=ff_filename,
                          **{i: j for i, j in kwargs.items() if i in
                             ilftr.RunAntechamber.required_params +
                             ilftr.RunAntechamber.optional_params +
                             ilftr.RunParmchk.required_params +
                             ilftr.RunParmchk.optional_params +
                             ilftwi.WriteTleapScript.required_params +
                             ilftwi.WriteTleapScript.optional_params +
                             ilftr.RunTleap.required_params +
                             ilftr.RunTleap.optional_params +
                             ilfpo.ProcessPrmtop.required_params +
                             ilfpo.ProcessPrmtop.optional_params})

        elif operation_type == "get_from_prmtop":
            if not isinstance(data, str):
                raise TypeError('data must be a str of the path to the prmtop '
                                'file for operation_type="get_from_prmtop"')
            tasks.append(ilfpo.ProcessPrmtop(
                            molecule=mol,
                            working_dir=working_dir,
                            db=db,
                            prmtop_path=data,
                            unique_molecule_name=label,
                            save_ff_to_db=save_ff_to_db,
                            save_ff_to_file=save_ff_to_file,
                            ff_filename=ff_filename,
                            **{i: j for i, j in kwargs.items() if i in
                               ilfpo.ProcessPrmtop.required_params
                               + ilfpo.ProcessPrmtop.optional_params}))

        elif operation_type == "get_from_dict":
            if not isinstance(data, dict):
                raise TypeError('data must be a dict for '
                                'operation_type="get_from_dict"')
            tasks.append(ilftwi.LabelFFDict(
                            mol=mol,
                            unlabeled_dict=data,
                            working_dir=working_dir,
                            label=label,
                            **{i: j for i, j in kwargs.items() if i in
                               ilftwi.LabelFFDict.required_params +
                               ilftwi.LabelFFDict.optional_params}))

        elif operation_type == "get_from_file":
            if not isinstance(data, str):
                raise TypeError('data must be a str of the path to the ff '
                                'json file for operation_type="get_from_file"')
            tasks.append(ilftwi.LabelFFDict(
                            ff_file=data,
                            working_dir=working_dir,
                            label=label,
                            **{i: j for i, j in kwargs.items() if i in
                               ilftwi.LabelFFDict.required_params +
                               ilftwi.LabelFFDict.optional_params}))

        elif operation_type == "get_from_db":
            tasks.append(ilftwi.LabelFFDictFromDB(
                            filter=data,
                            db=db,
                            working_dir=working_dir,
                            label=label,
                            **{i: j for i, j in kwargs.items()
                               if i in ilftwi.LabelFFDictFromDB.required_params
                               + ilftwi.LabelFFDictFromDB.optional_params}))

        # print(tasks)

        super(GetFFDictFW, self).__init__(tasks,
                                          parents=parents,
                                          name=name,
                                          **{i: j for i, j in kwargs.items()
                                             if i in FIREWORK_KWARGS})


class RunLammpsFW(fw.Firework):
    def __init__(self,
                 control_file=None,
                 db=None,
                 name="run_lammps",
                 parents=None,
                 working_dir=None,
                 save_run_to_db=True,
                 save_run_to_file=False,
                 **kwargs):
        """"""
        tasks = []
        working_dir = working_dir or os.getcwd()

        if not control_file:
            tasks.append(ilftwi.WriteControlFile(
                            working_dir=working_dir,
                            db=db,
                            save_to_db=save_run_to_db,
                            save_to_file=save_run_to_file,
                            **{i: j for i, j in kwargs.items() if i in
                               ilftwi.WriteControlFile.required_params
                               + ilftwi.WriteControlFile.optional_params}))

        tasks.append(ilftr.RunLammps(working_dir = working_dir,
                               **{i: j for i, j in kwargs.items() if i in
                                  ilftr.RunLammps.required_params + ilftr.RunLammps.optional_params}))

        super(RunLammpsFW, self).__init__(tasks,
                                          parents=parents,
                                          name=name,
                                          **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS})


if __name__ == "__main__":
    pass
