# coding: utf-8


# Defines standard lammps workflows.

import os

from fireworks.core.firework import Firework, Workflow

from mispr.gaussian.workflows.base.core import _process_mol_check
from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.lammps.defaults import (
    NPT_SETTINGS,
    NVT_SETTINGS,
    EMIN_SETTINGS,
    LAMMPS_RECIPE,
    RECIPE_SETTINGS,
    GENERAL_QADAPTER,
    QADAPTER_RUN_LAMMPS_SPEC,
)
from mispr.lammps.fireworks.core import GetFFDictFW, RunLammpsFW
from mispr.lammps.firetasks.write_inputs import WriteDataFile
from mispr.lammps.firetasks.parse_outputs import (
    GetRDF,
    CalcCN,
    CalcDiff,
    ExtractClusters,
    ProcessAnalysis,
)

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
    tag="unknown",
    **kwargs,
):
    """
    :param system_species_data: [dict]
              {
              species1_label: {"molecule": any format supported by mispr,
                               "molecule_operation_type": any format supported by mispr,
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
    :return:
    """
    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    system_mixture_data = {}

    if system_mixture_type == "concentration":
        system_mixture_data = {"Solutes": {}, "Solvents": {}}

    for label, ff_data in system_species_data.items():
        # Generate ff files in separate directories
        mol = recursive_relative_to_absolute_path(ff_data["molecule"], working_dir)
        mol_charge = ff_data.get("charge", None)
        mol_operation_type, mol, label, _, _, _ = _process_mol_check(
            working_dir=working_dir,
            mol_operation_type=ff_data["molecule_operation_type"],
            mol=mol,
            mol_name=label,
            db=db,
            charge=mol_charge,
            process_mol_func=kwargs.get("process_mol_func", False),
        )
        # TODO: add logic for ensuring that directory name is legal or at \
        #  least usable
        cur_species_dir = os.path.join(working_dir, label)

        # TODO: add logic for not setting file names that do not exist, \
        #  but code should work as is
        mol2_filename = f"{label}.mol2"
        frcmod_filename = f"{label}.frcmod"
        prmtop_filename = f"{label}.prmtop"
        inpcrd_filename = f"{label}.inpcrd"

        save_ff_to_db = ff_data.get("save_ff_to_db")
        save_ff_to_file = ff_data.get("save_ff_to_file")

        ff_fw = GetFFDictFW(
            mol=mol,
            mol_operation_type=mol_operation_type,
            data=ff_data["ff_param_data"],
            operation_type=ff_data["ff_param_method"],
            label=label,
            db=db,
            working_dir=cur_species_dir,
            output_filename_a=mol2_filename,
            input_filename_p=mol2_filename,
            output_filename_p=frcmod_filename,
            prmtop_filename=prmtop_filename,
            save_ff_to_db=save_ff_to_db,
            save_ff_to_file=save_ff_to_file,
            tleap_settings={
                "mol2_file_path": os.path.join(cur_species_dir, mol2_filename),
                "frcmod_file_path": os.path.join(cur_species_dir, frcmod_filename),
                "prmtop_file_path": os.path.join(cur_species_dir, prmtop_filename),
                "inpcrd_file_path": os.path.join(cur_species_dir, inpcrd_filename),
            },
            charge=mol_charge,
            **kwargs,
        )
        fireworks.append(ff_fw)

        if system_mixture_type == "concentration":
            system_mixture_data[ff_data["mol_mixture_type"]][label] = ff_data[
                "mixture_data"
            ]
        elif system_mixture_type == "number of molecules":
            system_mixture_data[label] = ff_data["mixture_data"]

    spec = kwargs.pop("spec", {})
    spec.update({"tag": tag, "_launch_dir": working_dir})
    extra_data_file_name = kwargs.pop("data_filename", "")
    data_fw = Firework(
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
        name="write_data_file",
        parents=fireworks[-1],
        spec=spec,
    )
    fireworks.append(data_fw)
    return fireworks


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

    default_recipe_step_names = [step[0] for step in LAMMPS_RECIPE]

    for index, step in enumerate(recipe):
        cur_step_dir = os.path.join(working_dir, step[0])
        lammpsin_key = f"run_{index}"
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
                lammpsin_key=lammpsin_key,
                **kwargs,
            )
            fireworks.append(cur_firework)
        elif step[1][0] == "template_str":
            cur_setting = recipe_settings[index]
            cur_firework = RunLammpsFW(
                working_dir=cur_step_dir,
                template_str=step[1][1],
                control_settings=cur_setting,
                db=db,
                lammpsin_key=lammpsin_key,
                spec=cur_spec_dict,
                **kwargs,
            )
            fireworks.append(cur_firework)
    return fireworks


def lammps_analysis_fws(analysis_list, analysis_settings, working_dir, **kwargs):
    fireworks = []
    links_dict = {}
    analysis_fw_ids = {}
    analysis_dir = os.path.join(working_dir, "analysis")

    if len(analysis_list) != len(analysis_settings):
        raise ValueError(
            f"{len(analysis_list)} types of analysis are requested while "
            f"{len(analysis_settings)} types of analysis settings are "
            f"provided"
        )

    for index, analysis in enumerate(analysis_list):
        cur_settings = analysis_settings[index]
        if analysis == "diffusion":
            diff_dir = os.path.join(analysis_dir, "diff")
            name = "diffusion_analysis"
            cur_settings.update(
                {
                    "outputs_dir": os.path.abspath(
                        os.path.join(diff_dir, "..", "..", "nvt")
                    )
                }
            )
            cur_firework = Firework(
                CalcDiff(diff_settings=cur_settings, working_dir=diff_dir),
                name=name,
                spec={"_launch_dir": diff_dir},
            )
            fireworks.append(cur_firework)

        elif analysis == "rdf":
            rdf_dir = os.path.join(analysis_dir, "rdf")
            name = "rdf_analysis"
            cur_settings.update(
                {
                    "filename": os.path.abspath(
                        os.path.join(rdf_dir, "..", "..", "nvt", "dump.nvt.*.dump")
                    )
                }
            )
            cur_firework = Firework(
                GetRDF(rdf_settings=cur_settings, working_dir=rdf_dir),
                name=name,
                spec={"_launch_dir": rdf_dir},
            )
            fireworks.append(cur_firework)

        elif analysis == "cn":
            cn_dir = os.path.join(analysis_dir, "cn")
            name = "cn_analysis"
            cur_settings.update(
                {
                    "filename": os.path.abspath(
                        os.path.join(cn_dir, "..", "..", "nvt", "dump.nvt.*.dump")
                    )
                }
            )
            cur_firework = Firework(
                CalcCN(cn_settings=cur_settings, working_dir=cn_dir),
                name=name,
                spec={"_launch_dir": cn_dir},
            )
            fireworks.append(cur_firework)

        elif analysis == "clusters":
            clusters_dir = os.path.join(analysis_dir, "clusters")
            name = "cluster_analysis"
            cur_settings.update(
                {
                    "filename": os.path.abspath(
                        os.path.join(clusters_dir, "..", "..", "nvt", "dump.nvt.*.dump")
                    )
                }
            )
            cur_firework = Firework(
                ExtractClusters(
                    cluster_settings=cur_settings, working_dir=clusters_dir
                ),
                name=name,
                spec={"_launch_dir": clusters_dir},
            )
            fireworks.append(cur_firework)

        else:
            raise ValueError(f"Unsupported analysis type: {analysis}")

        analysis_fw_ids[name] = [cur_firework.fw_id, cur_settings]

    if "cn_analysis" in analysis_fw_ids:
        if "r_cut" not in analysis_fw_ids["cn_analysis"][1]:
            if "rdf_analysis" not in analysis_fw_ids:
                raise ValueError(
                    "Cutoff distance required for calculating the CN is "
                    "not found and cannot be computed without performing "
                    "an RDF analysis"
                )
            else:
                links_dict[analysis_fw_ids["rdf_analysis"][0]] = analysis_fw_ids[
                    "cn_analysis"
                ][0]
        if "cluster_analysis" in analysis_fw_ids:
            if "r_cut" not in analysis_fw_ids["cluster_analysis"][1]:
                links_dict[analysis_fw_ids["cn_analysis"][0]] = analysis_fw_ids[
                    "cluster_analysis"
                ][0]

    if "cluster_analysis" in analysis_fw_ids and "cn_analysis" not in analysis_fw_ids:
        if "r_cut" not in analysis_fw_ids["cluster_analysis"][1]:
            raise ValueError(
                "Cutoff distance required for extracting clusters is "
                "not found and cannot be computed without performing "
                "a CN analysis"
            )

    final_analysis_fw = Firework(
        ProcessAnalysis(
            analysis_list=analysis_list,
            working_dir=analysis_dir,
            **{
                i: j
                for i, j in kwargs.items()
                if i
                in ProcessAnalysis.required_params + ProcessAnalysis.optional_params
            },
        ),
        parents=fireworks[:],
        name="analysis_postprocessing",
        spec={"_launch_dir": analysis_dir},
    )
    fireworks.append(final_analysis_fw)
    return fireworks, links_dict


def lammps_workflow(
    system_species_data=None,
    system_mixture_type=None,
    box_data=None,
    box_data_type="cubic",
    data_file_name="data.mixture",
    recipe=LAMMPS_RECIPE,
    recipe_settings=RECIPE_SETTINGS,
    recipe_qadapter=QADAPTER_RUN_LAMMPS_SPEC,
    db=None,
    working_dir=None,
    analysis_list=None,
    analysis_settings=None,
    name="lammps_workflow",
    **kwargs,
):

    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    links_dict = {}
    if all([system_species_data, system_mixture_type, box_data]):
        data_fws = lammps_data_fws(
            system_species_data,
            system_mixture_type,
            box_data,
            box_data_type=box_data_type,
            data_file_name=data_file_name,
            working_dir=working_dir,
            db=db,
            **kwargs,
        )
        fireworks += data_fws

    run_fws = lammps_run_fws(
        recipe=recipe,
        recipe_settings=recipe_settings,
        recipe_qadapter=recipe_qadapter,
        db=db,
        working_dir=working_dir,
        **kwargs,
    )
    fireworks += run_fws

    for index, fw in enumerate(fireworks[:-1]):
        links_dict[fw.fw_id] = fireworks[index + 1].fw_id

    if analysis_list:
        if not analysis_settings:
            analysis_settings = [{}] * len(analysis_list)
        analysis_fws, analysis_links_dict = lammps_analysis_fws(
            analysis_list, analysis_settings, working_dir=working_dir, **kwargs
        )
        links_dict[fireworks[-1].fw_id] = [i.fw_id for i in analysis_fws]
        fireworks += analysis_fws
        links_dict.update(analysis_links_dict)

    return Workflow(fireworks, links_dict=links_dict, name=name)
