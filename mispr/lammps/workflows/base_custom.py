# coding: utf-8


# Defines custom lammps workflows.

import os

from fireworks.core.firework import Firework, Workflow

from mispr.lammps.defaults import (
    NPT_SETTINGS,
    NVT_SETTINGS,
    EMIN_SETTINGS,
    LAMMPS_RECIPE,
    RECIPE_SETTINGS,
    GENERAL_QADAPTER,
    QADAPTER_RUN_LAMMPS_SPEC,
)
from mispr.lammps.fireworks.core_custom import GetFFDictFW, RunLammpsFW
from mispr.lammps.firetasks.write_inputs import WriteDataFile

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"


def write_lammps_data(
    system_species_data,
    system_mixture_type,
    box_data,
    box_data_type="cubic",
    data_file_name="complex.data",
    working_dir=None,
    db=None,
    **kwargs
):
    """
    :param system_species_data: [dict]
              {
              species1_label: {"molecule": pmg.Molecule,
                               "ff_param_method": str ("get_from_esp",
                                                       "get_from_prmtop",
                                                       "get_from_dict";
                                                       defaults to "get_from_esp"),
                               "ff_param_data": str (path_to_file) or dict (unlabeled_ff_dict),
                               "mol_mixture_type": "Solutes" or "Solvents",
                               "mixture_data": int or dict (depends on system_mixture_type)},
              ...,
              speciesN_label: {...}
              }
    :param system_mixture_type: [str] "concentration" or "number of molecules"
    :param box_data: [float, int, list (3,2), array (3,2), or LammpsBox]
        Definitions for box size. See box_data_type for info how to
        define this parameter.
    :param box_data_type: [string] Can be one of the following: 'cubic',
        'rectangular', or 'LammpsBox'. If 'cubic', box_data must be a
        float or int; if 'rectangular', box_data must be an array-like
        with size (3,2); if 'LammpsBox', box_data must be a
        pymatgen.io.lammps.data.LammpsBox object. Defaults to 'cubic'
    :param data_file_name:
    :param working_dir:
    :param db:
    :return:
    """
    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    if system_mixture_type == "concentration":
        system_mixture_data = {"Solutes": {}, "Solvents": {}}
    elif system_mixture_type == "number of molecules":
        system_mixture_data = {}

    links_dict = {}

    for species, ff_data in system_species_data.items():
        # Generate ff files in separate directories
        # TODO: add logic for ensuring that directory name is legal or at least usable
        cur_species_dir = os.path.join(working_dir, species)
        os.makedirs(cur_species_dir, exist_ok=True)

        # TODO: add logic for not setting file names that do not exist, but code should
        #  work as is
        mol2_filename = "{}.mol2".format(species)
        frcmod_filename = "{}.frcmod".format(species)
        prmtop_filename = "{}.prmtop".format(species)
        inpcrd_filename = "{}.inpcrd".format(species)

        cur_firework = GetFFDictFW(
            ff_data["molecule"],
            ff_data["ff_param_data"],
            operation_type=ff_data["ff_param_method"],
            label=species,
            db=db,
            working_dir=cur_species_dir,
            output_filename_a=mol2_filename,
            input_filename_p=mol2_filename,
            output_filename_p=frcmod_filename,
            prmtop_filename=prmtop_filename,
            tleap_settings={
                "mol2_file_path": os.path.join(cur_species_dir, mol2_filename),
                "frcmod_file_path": os.path.join(cur_species_dir, frcmod_filename),
                "prmtop_file_path": os.path.join(cur_species_dir, prmtop_filename),
                "inpcrd_file_path": os.path.join(cur_species_dir, inpcrd_filename),
            },
            **kwargs
        )

        fireworks.append(cur_firework)

        links_dict[cur_firework.fw_id] = []
        if system_mixture_type == "concentration":
            system_mixture_data[ff_data["mol_mixture_type"]][species] = ff_data[
                "mixture_data"
            ]
        elif system_mixture_type == "number of molecules":
            system_mixture_data[species] = ff_data["mixture_data"]

    spec = kwargs.pop("spec", {})
    firework2 = Firework(
        [
            WriteDataFile(
                working_dir=working_dir,
                data_filename=data_file_name,
                system_mixture_data_type=system_mixture_type,
                system_mixture_data=system_mixture_data,
                system_box_data=box_data,
                system_box_data_type=box_data_type,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in WriteDataFile.required_params + WriteDataFile.optional_params
                }
            )
        ],
        name="Write_Data_File_FW",
        parents=fireworks[-1],
        spec=spec,
    )

    fireworks.append(firework2)
    firework_ids = list(links_dict.keys())
    for index, firework_id in enumerate(firework_ids[:-1]):
        links_dict[firework_id].append(firework_ids[index + 1])
    links_dict[firework_ids[-1]].append(firework2.fw_id)
    return Workflow(fireworks, links_dict)


def run_lammps_recipe(
    recipe=LAMMPS_RECIPE,
    recipe_settings=RECIPE_SETTINGS,
    recipe_qadapter=QADAPTER_RUN_LAMMPS_SPEC,
    working_dir=None,
    db=None,
    init_spec=None,
    **kwargs
):
    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    links_dict = {}

    for index, step in enumerate(recipe):
        cur_step_dir = os.path.join(working_dir, step[0])
        os.makedirs(cur_step_dir, exist_ok=True)

        if index < len(QADAPTER_RUN_LAMMPS_SPEC):
            cur_qadapter_spec = QADAPTER_RUN_LAMMPS_SPEC[index].copy()
        else:
            cur_qadapter_spec = GENERAL_QADAPTER
        cur_qadapter_spec.update(recipe_qadapter[index])

        if step[1][0] == "template_filename":
            if step[1][1] == "emin_gaff":
                cur_setting = EMIN_SETTINGS.copy()
            elif step[1][1] == "npt":
                cur_setting = NPT_SETTINGS.copy()
            elif step[1][1] == "nvt":
                cur_setting = NVT_SETTINGS.copy()
            else:
                cur_setting = recipe_settings[index]
            cur_setting.update(recipe_settings[index])
            cur_spec_dict = {}
            if index == 0 and init_spec is not None:
                cur_spec_dict = init_spec.copy()
            cur_spec_dict.update({"_queueadapter": cur_qadapter_spec})
            cur_firework = RunLammpsFW(
                working_dir=cur_step_dir,
                template_filename=step[1][1],
                control_settings=cur_setting,
                db=db,
                spec=cur_spec_dict,
                **kwargs
            )
        elif step[1][0] == "template_str":
            cur_setting = recipe_settings[index]
            cur_firework = RunLammpsFW(
                working_dir=cur_step_dir,
                template_str=step[1][1],
                control_settings=cur_setting,
                db=db,
                spec={"_queueadapter": cur_qadapter_spec},
                **kwargs
            )

        fireworks.append(cur_firework)
        links_dict[cur_firework.fw_id] = []

    firework_ids = list(links_dict.keys())

    for index, firework_id in enumerate(firework_ids[:-1]):
        links_dict[firework_id].append(firework_ids[index + 1])

    return Workflow(fireworks, links_dict)


def fluid_workflow(
    system_species_data,
    system_mixture_type,
    box_data,
    box_data_type="cubic",
    data_file_name="complex.data",
    recipe=LAMMPS_RECIPE,
    recipe_settings=RECIPE_SETTINGS,
    recipe_qadapter=QADAPTER_RUN_LAMMPS_SPEC,
    db=None,
    working_dir=None,
    **kwargs
):
    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    if system_mixture_type == "concentration":
        system_mixture_data = {"Solutes": {}, "Solvents": {}}
    elif system_mixture_type == "number of molecules":
        system_mixture_data = {}

    links_dict = {}

    for species, ff_data in system_species_data.items():
        # Generate ff files in separate directories
        # TODO: add logic for ensuring that directory name is legal or at least usable
        cur_species_dir = os.path.join(working_dir, species)
        os.makedirs(cur_species_dir, exist_ok=True)

        # TODO: add logic for not setting file names that do not exist, but code should
        #  work as is
        mol2_filename = "{}.mol2".format(species)
        frcmod_filename = "{}.frcmod".format(species)
        prmtop_filename = "{}.prmtop".format(species)
        inpcrd_filename = "{}.inpcrd".format(species)

        cur_firework = GetFFDictFW(
            ff_data["molecule"],
            ff_data["ff_param_data"],
            operation_type=ff_data["ff_param_method"],
            label=species,
            db=db,
            working_dir=cur_species_dir,
            output_filename_a=mol2_filename,
            input_filename_p=mol2_filename,
            output_filename_p=frcmod_filename,
            prmtop_filename=prmtop_filename,
            tleap_settings={
                "mol2_file_path": os.path.join(cur_species_dir, mol2_filename),
                "frcmod_file_path": os.path.join(cur_species_dir, frcmod_filename),
                "prmtop_file_path": os.path.join(cur_species_dir, prmtop_filename),
                "inpcrd_file_path": os.path.join(cur_species_dir, inpcrd_filename),
            },
            **kwargs
        )

        fireworks.append(cur_firework)

        links_dict[cur_firework.fw_id] = []
        if system_mixture_type == "concentration":
            system_mixture_data[ff_data["mol_mixture_type"]][species] = ff_data[
                "mixture_data"
            ]
        elif system_mixture_type == "number of molecules":
            system_mixture_data[species] = ff_data["mixture_data"]

    spec = kwargs.pop("spec", {})
    firework2 = Firework(
        [
            WriteDataFile(
                working_dir=working_dir,
                data_filename=data_file_name,
                system_mixture_data_type=system_mixture_type,
                system_mixture_data=system_mixture_data,
                system_box_data=box_data,
                system_box_data_type=box_data_type,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in WriteDataFile.required_params + WriteDataFile.optional_params
                }
            )
        ],
        name="Write_Data_File_FW",
        parents=fireworks[-1],
        spec=spec,
    )

    fireworks.append(firework2)
    links_dict[firework2.fw_id] = []
    if not working_dir:
        working_dir = os.getcwd()
    run_fireworks = []

    for index, step in enumerate(recipe):
        cur_step_dir = os.path.join(working_dir, step[0])
        os.makedirs(cur_step_dir, exist_ok=True)
        if index < len(QADAPTER_RUN_LAMMPS_SPEC):
            cur_qadapter_spec = QADAPTER_RUN_LAMMPS_SPEC[index].copy()
        else:
            cur_qadapter_spec = GENERAL_QADAPTER.copy()
        cur_qadapter_spec.update(recipe_qadapter[index])
        if step[1][0] == "template_filename":
            if step[1][1] == "emin_gaff":
                cur_setting = EMIN_SETTINGS.copy()
            elif step[1][1] == "npt":
                cur_setting = NPT_SETTINGS.copy()
            elif step[1][1] == "nvt":
                cur_setting = NVT_SETTINGS.copy()
            else:
                cur_setting = recipe_settings[index]
            cur_setting.update(recipe_settings[index])
            cur_firework = RunLammpsFW(
                working_dir=cur_step_dir,
                db=db,
                template_filename=step[1][1],
                control_settings=cur_setting,
                spec={"_queueadapter": cur_qadapter_spec},
                **kwargs
            )
        elif step[1][0] == "template_str":
            cur_setting = recipe_settings[index]
            cur_firework = RunLammpsFW(
                working_dir=cur_step_dir,
                db=db,
                template_str=step[1][1],
                control_settings=cur_setting,
                spec={"_queueadapter": cur_qadapter_spec},
                **kwargs
            )

        run_fireworks.append(cur_firework)
        links_dict[cur_firework.fw_id] = []

    fireworks += run_fireworks
    firework_ids = list(links_dict.keys())

    for index, firework_id in enumerate(firework_ids[:-1]):
        links_dict[firework_id].append(firework_ids[index + 1])
    return Workflow(fireworks, links_dict)


if __name__ == "__main__":
    emin_settings = {
        "data_file_name": "../complex.data",
        "restart_final_name": "restart.emin.restart",
    }
    print(emin_settings)
    print(EMIN_SETTINGS)
    new_settings = EMIN_SETTINGS.copy()
    new_settings.update(emin_settings)
    print(new_settings)
    print(EMIN_SETTINGS)
