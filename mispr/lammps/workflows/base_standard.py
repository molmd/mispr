# coding: utf-8


# Defines standard lammps workflows.

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
from mispr.lammps.firetasks.write_inputs import WriteDataFile
from mispr.lammps.firetasks.parse_outputs import GetMSD, GetRDF, CalcDiff
from mispr.lammps.fireworks.core_standard import GetFFDictFW, RunLammpsFW

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Apr 2020"
__version__ = "0.0.1"


def lammps_data_fws(
    system_species_data,
    system_mixture_type,
    box_data,
    box_data_type="cubic",
    data_file_name="complex.data",
    working_dir=None,
    db=None,
    order_fws=False,
    **kwargs,
):
    """
    :param system_species_data: [dict]
              {
              species1_label: {"molecule": pmg.Molecule,
                               "ff_param_method": str can be one of the
                                    following ("get_from_esp",
                                               "get_from_prmtop",
                                               "get_from_dict")
                                    defaults to "get_from_esp",
                               "ff_param_data": str (path_to_file)
                                                or dict (unlabeled_ff_dict),
                               "mol_mixture_type": "Solutes" or "Solvents",
                               "mixture_data": int or dict depends on
                                    system_mixture_type},
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
    :param order_fws:
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
        # TODO: add logic for ensuring that directory name is legal or at \
        #  least usable
        cur_species_dir = os.path.join(working_dir, species)

        # TODO: add logic for not setting file names that do not exist, \
        #  but code should work as is
        mol2_filename = f"{species}.mol2"
        frcmod_filename = f"{species}.frcmod"
        prmtop_filename = f"{species}.prmtop"
        inpcrd_filename = f"{species}.inpcrd"

        # TODO: Add ability to not need to supply db and still save to db.
        #  Change defaults to those under "if db:" only.
        if db:
            save_ff_to_db = ff_data.get("save_ff_to_db", True)
            save_ff_to_file = ff_data.get("save_ff_to_file", False)
        else:
            save_ff_to_db = ff_data.get("save_ff_to_db", False)
            save_ff_to_file = ff_data.get("save_ff_to_file", True)

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
            **kwargs,
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
    extra_data_file_name = kwargs.pop("data_filename", "")
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
                },
            )
        ],
        name="Write_Data_File_FW",
        parents=fireworks[-1],
        spec=spec,
    )

    fireworks.append(firework2)
    firework_ids = list(links_dict.keys())
    if order_fws:
        for index, firework_id in enumerate(firework_ids[:-1]):
            links_dict[firework_id].append(firework_ids[index + 1])
        links_dict[firework_ids[-1]].append(firework2.fw_id)

    links_dict[firework2.fw_id] = []
    return fireworks, links_dict


def lammps_run_fws(
    recipe=LAMMPS_RECIPE,
    recipe_settings=RECIPE_SETTINGS,
    recipe_qadapter=QADAPTER_RUN_LAMMPS_SPEC,
    init_spec=None,
    db=None,
    working_dir=None,
    save_runs_to_db=True,
    save_runs_to_file=False,
    **kwargs,
):

    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    links_dict = {}

    default_recipe_step_names = [step[0] for step in LAMMPS_RECIPE]

    for index, step in enumerate(recipe):
        cur_step_dir = os.path.join(working_dir, step[0])

        matching_default_step_index = None
        if step[0] in default_recipe_step_names:
            matching_default_step_index = default_recipe_step_names.index(step[0])
            cur_qadapter_spec = QADAPTER_RUN_LAMMPS_SPEC[
                matching_default_step_index
            ].copy()
        else:
            cur_qadapter_spec = GENERAL_QADAPTER.copy()
        cur_qadapter_spec.update(recipe_qadapter[index])

        cur_spec_dict = {}
        if index == 0 and init_spec is not None:
            cur_spec_dict = init_spec.copy()
        cur_spec_dict.update({"_queueadapter": cur_qadapter_spec})

        if step[1][0] == "template_filename":
            if step[1][1] == "emin_gaff" or step[1][1] == "emin":
                cur_setting = EMIN_SETTINGS.copy()
            elif step[1][1] == "npt":
                cur_setting = NPT_SETTINGS.copy()
            elif step[1][1] == "nvt":
                cur_setting = NVT_SETTINGS.copy()
            else:
                cur_setting = recipe_settings[index].copy()
            cur_setting.update(recipe_settings[index])

            cur_firework = RunLammpsFW(
                working_dir=cur_step_dir,
                db=db,
                template_filename=step[1][1],
                control_settings=cur_setting,
                spec=cur_spec_dict,
                save_run_to_db=save_runs_to_db,
                save_run_to_file=save_runs_to_file,
                **kwargs,
            )
        elif step[1][0] == "template_str":
            cur_setting = recipe_settings[index]
            cur_firework = RunLammpsFW(
                working_dir=cur_step_dir,
                template_str=step[1][1],
                control_settings=cur_setting,
                db=db,
                spec=cur_spec_dict,
                **kwargs,
            )
        fireworks.append(cur_firework)
        links_dict[cur_firework.fw_id] = []
    return fireworks, links_dict


def lammps_analysis_fws(analysis_list, analysis_settings, working_dir=None, **kwargs):
    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    links_dict = {}

    for index, type in enumerate(analysis_list):
        if type == "msd_from_dump":
            cur_firework = Firework(
                GetMSD(
                    msd_method="from_dump",
                    msd_settings=analysis_settings[index],
                    working_dir=os.path.join(working_dir, "msd"),
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i in GetMSD.required_params + GetMSD.optional_params
                    },
                )
            )
            fireworks.append(cur_firework)
            links_dict[cur_firework.fw_id] = []

        elif type == "msd_from_log":
            cur_firework = Firework(
                GetMSD(
                    msd_method="from_log",
                    msd_settings=analysis_settings[index],
                    working_dir=os.path.join(working_dir, "msd"),
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i in GetMSD.required_params + GetMSD.optional_params
                    },
                )
            )
            fireworks.append(cur_firework)
            links_dict[cur_firework.fw_id] = []

        elif type == "diffusion":
            cur_firework = Firework(
                CalcDiff(
                    diff_settings=analysis_settings[index],
                    working_dir=os.path.join(working_dir, "diff"),
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i in CalcDiff.required_params + CalcDiff.optional_params
                    },
                )
            )
            fireworks.append(cur_firework)
            links_dict[cur_firework.fw_id] = []

        elif type == "rdf":
            cur_working_dir = os.path.join(working_dir, "rdf")

            cur_settings = analysis_settings[index].copy()
            cur_settings.update(
                {
                    "filename": os.path.abspath(
                        os.path.join(cur_working_dir, "..", "nvt", "dump.nvt.*.dump")
                    )
                }
            )
            cur_firework = Firework(
                GetRDF(
                    rdf_settings=cur_settings,
                    working_dir=cur_working_dir,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i in GetRDF.required_params + GetRDF.optional_params
                    },
                )
            )
            fireworks.append(cur_firework)
            links_dict[cur_firework.fw_id] = []

        else:
            # TODO: raise exception
            pass

    return fireworks, links_dict


def lammps_workflow(
    system_species_data=None,
    system_mixture_type=None,
    box_data=None,
    box_data_type="cubic",
    data_file_name="complex.data",
    recipe=LAMMPS_RECIPE,
    recipe_settings=RECIPE_SETTINGS,
    recipe_qadapter=QADAPTER_RUN_LAMMPS_SPEC,
    db=None,
    working_dir=None,
    analysis_list=None,
    analysis_settings=None,
    **kwargs,
):

    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    links_dict = {}
    if all([system_species_data, system_mixture_type, box_data]):
        fireworks, links_dict = lammps_data_fws(
            system_species_data,
            system_mixture_type,
            box_data,
            box_data_type=box_data_type,
            data_file_name=data_file_name,
            working_dir=working_dir,
            db=db,
            **kwargs,
        )

    run_fireworks, run_links_dict = lammps_run_fws(
        recipe=recipe,
        recipe_settings=recipe_settings,
        recipe_qadapter=recipe_qadapter,
        db=db,
        working_dir=working_dir,
        **kwargs,
    )

    fireworks += run_fireworks
    links_dict.update(run_links_dict)

    if analysis_list:
        if not analysis_settings:
            analysis_settings = [{}] * len(analysis_list)
        analysis_fireworks, analysis_links_dict = lammps_analysis_fws(
            analysis_list, analysis_settings, working_dir=working_dir, **kwargs
        )
        fireworks += analysis_fireworks
        links_dict.update(analysis_links_dict)

    # After all fws set
    fw_ids = list(links_dict.keys())
    for index, fw_id in enumerate(fw_ids[:-1]):
        links_dict[fw_id].append(fw_ids[index + 1])
    return Workflow(fireworks, links_dict)


if __name__ == "__main__":
    import fireworks as fw

    launchpad = fw.LaunchPad(
        host="mongodb://superuser:<password>@localhost:27017/fireworks?authSource=admin",
        uri_mode=True,
    )
    launchpad.reset("", require_password=False)
