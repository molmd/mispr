"""Define standard lammps workflows."""

import os
import numpy as np

from fireworks.core.firework import Firework, Workflow

from mispr.lammps.defaults import (
    MIN_INTERFACE_SETTINGS,
    NVE_INTERFACE_SETTINGS,
    NVT_INTERFACE_SETTINGS,
    INTERFACE_LAMMPS_RECIPE,
    INTERFACE_RECIPE_SETTINGS,
    GENERAL_QADAPTER,
    INTERFACE_QADAPTER_RUN_LAMMPS_SPEC,
)
from mispr.lammps.fireworks.core import GetFFDictFW, RunLammpsFW
from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.workflows.base.core import _process_mol_check
from mispr.lammps.firetasks.write_inputs import (
    WriteDataFile,
    RunPackmol,
    AddSlabToXYZ,
)
from mispr.lammps.firetasks.parse_outputs import (
    CalcCN,
    GetRDF,
    CalcDiff,
    ExtractClusters,
    ProcessAnalysis,
)

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@stonybrook.edu"
__status__ = "Development"
__date__ = "Jun 2024"
__version__ = "0.0.4"


def lammps_interface_data_fws(
    system_species_data,
    system_mixture_data,
    box_data,
    system_mixture_type="number of molecules",
    box_data_type="rectangular",
    data_file_name="complex.data",
    bulk_molecule_names=None,
    bulk_molecule_files=None,
    bulk_num_molecules=None,
    bulk_box_data=None,
    slab_file=None,
    sandwich=True,
    working_dir=None,
    molecule_dir=None,
    db=None,
    tag="unknown",
    **kwargs,
):
    """
    Generate FireWorks for writing LAMMPS data files.

    Args:
        system_species_data (dict): Dictionary containing species data. The keys are the
            species labels and the values are dictionaries containing the following keys:

            * molecule (Molecule, GaussianOutput, str, dict): Source of the molecule to
              be processed. See ``process_mol`` in ``mispr/gaussian/utilities/mol.py``
              for supported operations.
            * molecule_operation_type (str): Type of molecule operation. Must match with
              ``molecule`` value.
            * ff_param_method (str): Method for obtaining force field parameters. Must
              match with ``ff_param_data`` value. See available methods for ``GetFFDictFw``
              in ``mispr/lammps/fireworks/core.py``.
            * ff_param_data (str or dict): Data regarding necessary information to obtain
              ff parameters. Must match with ``ff_param_method`` value.
            * mol_mixture_type (str): Type of mixture data. Must be "Solutes" or "Solvents".
            * mixture_data (int or dict): Information regarding the number of molecules
              of this type in the system. Depends on the ``system_mixture_type`` parameter.

        system_mixture_data (dict): Dictionary containing information related to the number
            of molecules of each species in the system. Currently the only supported type
            is a dictionary where the keys are the species labels and the values are the
            number of molecules of that type in the system. This corresponds with the
            "number of molecules" value for the ``system_mixture_type`` parameter.
            Defaults to "number of molecules".
        box_data (float, int, list (3,2), array (3,2), or LammpsBox): Definitions for
            box size. See ``box_data_type`` for info how to define this parameter.
        system_mixture_type (str): Type of mixture data. Must be "concentration" or
            "number of molecules". See ``LammpsDataWrapper`` in
            ``pymatgen/io/lammps/data.py`` for more information.
        box_data_type (str, optional): Determines the value of the ``box_data`` parameter.
            Can be one of the following: 'cubic', 'rectangular', or 'LammpsBox'. If
            'cubic', ``box_data`` must be a float or int; if 'rectangular', ``box_data``
            must be an array-like with size (3,2); if 'LammpsBox', ``box_data`` must be
            a ``LammpsBox`` object. Defaults to 'cubic'.
        data_file_name (str, optional): Name of the data file to be written. Defaults to
            "complex.data".
        working_dir (str, optional): Directory where the data files will be written.
            Defaults to the current working directory.
        db (str or dict, optional): Database credentials. Could be a string with the path
            to the database file or a dictionary with the database credentials. If none
            is provided, attempts to read the configuration files. Only used when
            ``save_ff_to_db`` is ``True``.
        tag (str, optional): Tag for the Fireworks. Defaults to "unknown".
        kwargs (keyword arguments): Additional keyword arguments.

    Returns:
        fireworks (list): List of FireWorks for writing LAMMPS data files.
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
        ff_data["ff_param_data"] = recursive_relative_to_absolute_path(
            ff_data["ff_param_data"], working_dir
        )
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

    # Run packmol (optional)
    if all([bulk_molecule_names, bulk_molecule_files, bulk_num_molecules]):
        run_packmol = True
        packmol_dir = os.path.join(working_dir, "packmol")
        os.makedirs(packmol_dir, exist_ok=True)

        # convert bulk_box_data to format required by RunPackmol
        # ie, changes from [[xlo, xhi], ...] to [xlo, ylo, zlo, xhi, ...]
        revised_bulk_box_data = np.transpose(bulk_box_data).flatten().tolist()

        packmol_spec = kwargs.pop("spec", {})
        packmol_spec.update({"tag": tag, "_launch_dir": packmol_dir})

        pack_fw = Firework(
            [
                RunPackmol(
                    molecule_files=bulk_molecule_files,
                    molecule_dir=molecule_dir,
                    working_dir=packmol_dir,
                    output_file="bulk.xyz",
                    box_dim=revised_bulk_box_data,
                    num_mols=bulk_num_molecules,
                    mol_list=bulk_molecule_names,
                )
            ],
            name="run_packmol",
            parents=fireworks[-1],
            spec=packmol_spec,
        )
        fireworks.append(pack_fw)

    # Add slabs to bulk (optional)
    if slab_file:
        slab_dir = os.path.join(working_dir, "slab")
        os.makedirs(slab_dir, exist_ok=True)

        slab_spec = kwargs.pop("spec", {})
        slab_spec.update({"tag": tag, "_launch_dir": slab_dir})

        slab_fw = Firework(
            [
                AddSlabToXYZ(
                    slab_file=slab_file,
                    working_dir=slab_dir,
                    sandwich=sandwich,
                    **{
                        i: j
                        for i, j in kwargs.items()
                        if i
                        in AddSlabToXYZ.required_params + AddSlabToXYZ.optional_params
                    },
                )
            ],
            name="make_interface",
            parents=fireworks[-1],
            spec=slab_spec,
        )
        fireworks.append(slab_fw)

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


def lammps_interface_run_fws(
    recipe=INTERFACE_LAMMPS_RECIPE,
    recipe_settings=INTERFACE_RECIPE_SETTINGS,
    recipe_qadapter=INTERFACE_QADAPTER_RUN_LAMMPS_SPEC,
    init_spec=None,
    db=None,
    working_dir=None,
    save_runs_to_db=True,
    save_runs_to_file=False,
    **kwargs,
):
    """
    Generate FireWorks for running LAMMPS simulations.

    Args:
        recipe (list, optional): List of lists containing the name of the step and the
            template filename or string for the LAMMPS input file. Defaults to
            INTERFACE_LAMMPS_RECIPE.
        recipe_settings (list, optional): List of dictionaries containing the settings
            for each step in the recipe. Defaults to INTERFACE_RECIPE_SETTINGS.
        recipe_qadapter (list, optional): List of dictionaries containing the settings
            for the queue adapter for each step in the recipe. Defaults to
            INTERFACE_QADAPTER_RUN_LAMMPS_SPEC.
        init_spec (dict, optional): Initial spec for the FireWorks. Defaults to ``None``.
        db (str or dict, optional): Database credentials. Could be a string with the path
            to the database file or a dictionary with the database credentials. If none
            is provided, attempts to read the configuration files. Only used when
            ``save_runs_to_db`` is ``True``.
        working_dir (str, optional): Directory where the data files will be written.
            Defaults to the current working directory.
        save_runs_to_db (bool, optional): Whether to save the runs to the database.
            Defaults to ``True``.
        save_runs_to_file (bool, optional): Whether to save the runs to a file. Defaults
            to ``False``.
        kwargs (keyword arguments): Additional keyword arguments.

    Returns:
        fireworks (list): List of FireWorks for running LAMMPS simulations.
    """

    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []

    default_recipe_step_names = [step[0] for step in INTERFACE_LAMMPS_RECIPE]

    for index, step in enumerate(recipe):
        print(step)
        cur_step_dir = os.path.join(working_dir, step[0])
        lammpsin_key = f"run_{index}"
        matching_default_step_index = None
        if step[0] in default_recipe_step_names:
            matching_default_step_index = default_recipe_step_names.index(step[0])
            cur_qadapter_spec = INTERFACE_QADAPTER_RUN_LAMMPS_SPEC[
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
            if step[1][1] == "min_interface":
                cur_setting = MIN_INTERFACE_SETTINGS.copy()
            elif step[1][1] == "nve_interface":
                cur_setting = NVE_INTERFACE_SETTINGS.copy()
            elif step[1][1] == "nvt_interface":
                cur_setting = NVT_INTERFACE_SETTINGS.copy()
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

# TODO: Change lammps_analysis_fws to work for the interfacial systems
def lammps_interface_analysis_fws(analysis_list, analysis_settings, working_dir, **kwargs):
    """
    Generate FireWorks for running LAMMPS analysis.

    Args:
        analysis_list (list): List of analysis types to perform. Supported types are:
            'diffusion', 'rdf', 'cn', and 'clusters'.
        analysis_settings (list): List of dictionaries containing the settings for each
            analysis type.
        working_dir (str): Directory where the data files will be written.
        kwargs (keyword arguments): Additional keyword arguments.

    Returns:
        tuple:
            - fireworks (list): List of FireWorks objects for running LAMMPS analysis.
            - links_dict (dict): Dictionary containing the links between the FireWorks.
    """
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


def lammps_interfacial_workflow(
    system_species_data=None,
    system_mixture_data=None,
    box_data=None,
    system_mixture_type="number of molecules",
    box_data_type="rectangular",
    data_file_name="data.mixture",
    bulk_molecule_names=None,
    bulk_molecule_files=None,
    bulk_num_molecules=None,
    bulk_box_data=None,
    molecule_dir=None,
    slab_file=None,
    sandwich=True,
    recipe=INTERFACE_LAMMPS_RECIPE,
    recipe_settings=INTERFACE_RECIPE_SETTINGS,
    recipe_qadapter=INTERFACE_QADAPTER_RUN_LAMMPS_SPEC,
    db=None,
    working_dir=None,
    analysis_list=None,
    analysis_settings=None,
    name="lammps_workflow",
    **kwargs,
):
    """
    Create a LAMMPS workflow.

    Args:
        system_species_data (dict, optional): Dictionary containing species data. Refer
            to the ``lammps_data_fws`` function for more information. Defaults to
            ``None``. If not provided, the workflow will not create any FireWorks for
            writing LAMMPS data files.
        system_mixture_type (str): Type of mixture data. Must be "concentration" or
            "number of molecules". See ``LammpsDataWrapper`` in
            ``pymatgen/io/lammps/data.py`` for more information. Defaults to ``None``.
            If not provided, the workflow will not create any FireWorks for writing
            LAMMPS data files.
        box_data (float, int, list (3,2), array (3,2), or LammpsBox): Definitions for
            box size. See ``lammps_data_fws`` for info on how to define this parameter.
            Defaults to ``None``. If not provided, the workflow will not create any
            FireWorks for writing LAMMPS data files.
        box_data_type (str, optional): Determines the value of the ``box_data``
            parameter. Defaults to 'rectangular'.
        data_file_name (str, optional): Name of the data file to be written. Defaults
            to 'data.mixture'.
        bulk_molecule_names (list, optional): List of molecule labels
            for the bulk system. The order of this list determines the
            order the molecules are added by packmol. Only used when the
            bulk system is to be generated using the ``RunPackmol``
            firetask. Defaults to ``None``.
        bulk_molecule_files (list, optional): List of molecule files for
            the bulk system. Supports any filetype that is supported by
            the ``Molecule`` object. Only used when the bulk system is
            to be generated by the ``RunPackmol`` firetask. Defaults to
            ``None``.
        bulk_num_molecules (dict, optional): Dict where the keys are the
            molecule labels and the values are the number of molecules
            of that type to be added to the bulk system. Only used when
            the bulk system is to be generated by the ``RunPackmol``
            firetask. Defaults to ``None``.
        bulk_box_data (list, optional): List of box dimensions for the
            bulk part of the system. Intended to be used by the
            ``RunPackmol`` firetask. Defaults to ``None``.
        molecule_dir (str, optional): Directory where the molecule files
            are stored. Intended to be used by the ``RunPackmol``
            firetask. Defaults to ``None``.
        slab_file (str, optional): Path to the slab xyz file. Intended
            to be used by the ``AddSlabToXYZ`` firetask. Defaults to
            ``None``.
        sandwich (bool, optional): Whether to add a second slab to the
            system in the sandwich configuration. If ``False`` then will
            create an interfacial system with only one slab. Intented to
            be used by the ``AddSlabToXYZ`` firetask. Defaults to
            ``True``.
        recipe (list, optional): List of lists containing the name of the step and the
            template filename or string for the LAMMPS input file. Defaults to
            INTERFACE_LAMMPS_RECIPE.
        recipe_settings (list, optional): List of dictionaries containing the settings
            for each step in the recipe. Defaults to INTERFACE_RECIPE_SETTINGS.
        recipe_qadapter (list, optional): List of dictionaries containing the settings
            for the queue adapter for each step in the recipe. Defaults to
            INTERFACE_QADAPTER_RUN_LAMMPS_SPEC.
        db (str or dict, optional): Database credentials. Could be a string with the
            path to the database file or a dictionary with the database credentials.
            If none is provided, attempts to read the configuration files. Only used when
            ``save_runs_to_db`` is ``True``. Defaults to ``None``.
        working_dir (str, optional): Directory where the data files will be written.
            Defaults to the current working directory.
        analysis_list (list, optional): List of analysis types to perform. Supported
            types are: 'diffusion', 'rdf', 'cn', and 'clusters'. Defaults to ``None``.
            If not provided, the workflow will not create any FireWorks for running
            LAMMPS analysis.
        analysis_settings (list, optional): List of dictionaries containing the settings
            for each analysis type. Defaults to ``None``. If not provided, the workflow
            will not create any FireWorks for running LAMMPS analysis.
        name (str, optional): Name of the workflow. Defaults to 'lammps_workflow'.
        kwargs (keyword arguments): Additional keyword arguments.

    Returns:
        Workflow
    """

    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    links_dict = {}
    if all([system_species_data, system_mixture_type, box_data]):
        data_fws = lammps_interface_data_fws(
            system_species_data,
            system_mixture_data,
            box_data,
            system_mixture_type=system_mixture_type,
            box_data_type=box_data_type,
            data_file_name=data_file_name,
            working_dir=working_dir,
            db=db,
            bulk_molecule_names=bulk_molecule_names,
            bulk_molecule_files=bulk_molecule_files,
            bulk_num_molecules=bulk_num_molecules,
            bulk_box_data=bulk_box_data,
            slab_file=slab_file,
            sandwich=sandwich,
            molecule_dir=molecule_dir,
            **kwargs,
        )
        fireworks += data_fws

    run_fws = lammps_interface_run_fws(
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
        analysis_fws, analysis_links_dict = lammps_interface_analysis_fws(
            analysis_list, analysis_settings, working_dir=working_dir, **kwargs
        )
        links_dict[fireworks[-1].fw_id] = [i.fw_id for i in analysis_fws]
        fireworks += analysis_fws
        links_dict.update(analysis_links_dict)

    return Workflow(fireworks, links_dict=links_dict, name=name)
