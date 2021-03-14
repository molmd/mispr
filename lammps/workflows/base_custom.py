import os
from shutil import copy2
from fireworks.core.firework import Firework, Workflow
from infrastructure.gaussian.utils.utils import get_mol_formula
from infrastructure.lammps.firetasks.write_inputs import WriteDataFile
from infrastructure.lammps.fireworks.core_custom import GetFFDictFW, RunLammpsFW


EMIN_SETTING_DEFAULTS = {"neigh_args": "2.0 bin",
                         "neigh_modify_args": "delay 0 every 1 check yes",
                         "pair_style": "lj/cut/coul/long",
                         "pair_style_args": 10.0,
                         "kspace_style": "pppm",
                         "kspace_style_args": "1.0e-4",
                         "data_file_name": "complex.data",
                         "restart_final_name": "restart.emin.restart",
                         "dump_file_name": "dump.emin.dump",
                         "dump_modify_args": ""}

NPT_SETTING_DEFAULTS = {"neigh_args": "2.0 bin",
                        "neigh_modify_args": "delay 0 every 1 check yes",
                        "special_bonds_style": "amber",
                        "special_bonds_value": "",
                        "bond_style": "harmonic",
                        "bond_style_args": "",
                        "angle_style": "harmonic",
                        "dihedral_style": "harmonic",
                        "improper_style": "cvff",
                        "pair_style": "lj/cut/coul/long",
                        "pair_style_args": 10.0,
                        "kspace_style": "pppm",
                        "kspace_style_args": "1.0e-4",
                        "restart_filename": "restart.emin.restart",
                        "group_definitions": "",
                        "temperature_initial": 298.15,
                        "velocity_seed": 250,
                        "shake_logic": "#",
                        "shake_group": "",
                        "shake_topologies": "",
                        "temperature_final": 298.15,
                        "temp_damp": 100.0,
                        "pressure_type": "iso",
                        "pressure_initial": 1.0,
                        "pressure_final": 1.0,
                        "pres_damp": 1000.0,
                        "thermo": 1000,
                        "timestep": 1,
                        "dump_period": 50000,
                        "dump_file_name": "dump.npt.*.dump",
                        "dump_modify_logic": "#",
                        "dump_modify_args": "",
                        "restart_period": 1000000,
                        "restart_intermediate_file_name": "restart.npt.*.restart",
                        "run": 2000000,
                        "restart_final_file_name": "restart.npt.restart"}

NVT_SETTING_DEFAULTS = {"neigh_args": "2.0 bin",
                        "neigh_modify_args": "delay 0 every 1 check yes",
                        "special_bonds_style": "amber",
                        "special_bonds_value": "",
                        "bond_style": "harmonic",
                        "bond_style_args": "",
                        "angle_style": "harmonic",
                        "dihedral_style": "harmonic",
                        "improper_style": "cvff",
                        "pair_style": "lj/cut/coul/long",
                        "pair_style_args": 10.0,
                        "kspace_style": "pppm",
                        "kspace_style_args": "1.0e-4",
                        "restart_file_name": "restart.npt.restart",
                        "group_definitions": "",
                        "shake_logic": "#",
                        "shake_group": "",
                        "shake_topologies": "",
                        "temperature_initial": 298.15,
                        "temperature_final": 298.15,
                        "temp_damp": 100.0,
                        "thermo": 5,
                        "compute_definitions": "",
                        "thermo_style_compute": "",
                        "timestep": 1,
                        "dump_period": 50000,
                        "dump_file_name": "dump.nvt.*.dump",
                        "dump_modify_logic": "#",
                        "dump_modify_args": "",
                        "restart_period": 1000000,
                        "restart_intermediate_file_name": "restart.nvt.*.restart",
                        "run": 5000000,
                        "restart_final_file_name": "restart.nvt.restart"}

# Format of list in recipe is:
# [name for step in recipe to be used as the name of the directory,
#  [control template parameter type, control template data],
#  dictionary for settings]

LAMMPS_RECIPE_DEFAULT = [["emin", ["template_filename", "emin_gaff"]],
                         ["npt", ["template_filename", "npt"]],
                         ["melt", ["template_filename", "nvt"]],
                         ["quench", ["template_filename", "nvt"]],
                         ["nvt", ["template_filename", "nvt"]]]

RECIPE_SETTING_DEFAULT = [{"data_file_name": "../complex.data",
                           "restart_final_name": "restart.emin.restart"},
                          {"restart_filename": "../emin/restart.emin.restart",
                           "restart_final_file_name": "restart.npt.restart"},
                          {"restart_filename": "../npt/restart.npt.restart",
                           "temperature_initial": 500.0,
                           "temperature_final": 500.0,
                           "run": 1000000,
                           "restart_final_file_name": "restart.melt_500K.restart"},
                          {"restart_file_name": "../melt/restart.melt_500K.restart",
                           "temperature_initial": 500.0,
                           "temperature_final": 298.15,
                           "run": 2000000,
                           "restart_final_file_name": "restart.quench_298-15K.restart"},
                          {"restart_file_name": "../quench/restart.quench_298-15K.restart",
                           "temperature_initial": 298.15,
                           "temperature_final": 298.15,
                           "run": 5000000,
                           "restart_final_file_name": "restart.nvt_5-ns.restart"}]

QADAPTER_RUN_LAMMPS_SPEC_DEFAULTS = [{"walltime": "00:10:00",
                                      "job_name": "emin"},
                                     {"walltime": "08:00:00",
                                      "job_name": "npt"},
                                     {"walltime": "03:00:00",
                                      "job_name": "melt"},
                                     {"walltime": "06:00:00",
                                      "job_name": "quench"},
                                     {"walltime": "15:00:00",
                                      "job_name": "nvt"}]

def write_lammps_data(system_species_data,
                      system_mixture_type,
                      box_data,
                      box_data_type = 'cubic',
                      data_file_name = "complex.data",
                      working_dir=None,
                      **kwargs):
    """

    :param system_species_data: [dict]
              {
              species1_label: {"molecule": pmg.Molecule,
                               "ff_param_method": str ("get_from_esp",
                                                       "get_from_prmtop",
                                                       "get_from_dict";
                                                       defaults to "get_from_esp"),
                               "ff_param_data": str (path to file) or dict (unlabeled ff_dict),
                               "mol_mixture_type": "Solutes" or "Solvents",
                               "mixture_data": int or dict},
              ...,
              speciesN_label: {...}
              }
    :param system_mixture_type: [str] "concentration" or "number of molecules"
    :param box_data: [float, int, list (3,2), array (3,2), or LammpsBox] Definitions for box size.
        See box_data_type for info how to define this parameter.
    :param box_data_type: [string] Can be one of the following: 'cubic', 'rectangular', or 'LammpsBox'.
        If 'cubic', box_data must be a float or int; if 'rectangular', box_data must be an array-like
        with size (3,2); if 'LammpsBox', box_data must be a pymatgen.io.lammps.data.LammpsBox object.
        Defaults to 'cubic'
    :return:
    """
    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    if system_mixture_type == "concentration":
        system_mixture_data = {"Solutes": {},
                               "Solvents": {}}
    elif system_mixture_type == "number of molecules":
        system_mixture_data = {}

    links_dict = {}

    for species, ff_data in system_species_data.items():
        # Generate ff files in separate directories
        # TODO: add logic for ensuring that directory name is legal or at least usable
        cur_species_dir = os.path.join(working_dir, species)
        os.makedirs(cur_species_dir, exist_ok = True)


        # TODO: add logic for not setting file names that do not exist, but code should work as is
        mol2_filename = "{}.mol2".format(species)
        frcmod_filename = "{}.frcmod".format(species)
        prmtop_filename = "{}.prmtop".format(species)
        inpcrd_filename = "{}.inpcrd".format(species)

        cur_firework = GetFFDictFW(ff_data["molecule"],
                                   ff_data["ff_param_data"],
                                   operation_type = ff_data["ff_param_method"],
                                   label = species,
                                   working_dir = cur_species_dir,
                                   output_filename_a = mol2_filename,
                                   input_filename_p = mol2_filename,
                                   output_filename_p = frcmod_filename,
                                   prmtop_filename = prmtop_filename,
                                   tleap_settings = {"mol2_file_path": os.path.join(cur_species_dir, mol2_filename),
                                                     "frcmod_file_path": os.path.join(cur_species_dir, frcmod_filename),
                                                     "prmtop_file_path": os.path.join(cur_species_dir, prmtop_filename),
                                                     "inpcrd_file_path": os.path.join(cur_species_dir, inpcrd_filename)},
                                   **kwargs)

        fireworks.append(cur_firework)

        links_dict[cur_firework.fw_id] = []
        # print(cur_firework.fw_id)
        if system_mixture_type == "concentration":
            system_mixture_data[ff_data["mol_mixture_type"]][species] = ff_data["mixture_data"]
        elif system_mixture_type == "number of molecules":
            system_mixture_data[species] = ff_data["mixture_data"]

    spec = kwargs.pop("spec", {})
    firework2 = Firework([WriteDataFile(working_dir = working_dir,
                                        data_filename = data_file_name,
                                        system_mixture_data_type = system_mixture_type,
                                        system_mixture_data = system_mixture_data,
                                        system_box_data = box_data,
                                        system_box_data_type = box_data_type,
                                        **{i: j for i, j in kwargs.items() if i in
                                           WriteDataFile.required_params +
                                           WriteDataFile.optional_params})],
                         name = "Write_Data_File_FW",
                         parents = fireworks[-1],
                         spec = spec)

    fireworks.append(firework2)

    # print("FIREWORKS:")
    # for firework in fireworks:
    #     # print(firework.fw_id)
    #     # print(firework.fw_id % len(fireworks) + 1)
    #     print(firework.tasks)

    firework_ids = list(links_dict.keys())

    for index, firework_id in enumerate(firework_ids[:-1]):
        links_dict[firework_id].append(firework_ids[index + 1])
        # links_dict[firework_id] += firework_ids[index + 1:]
        # links_dict[firework_id].append(firework2.fw_id)
    links_dict[firework_ids[-1]].append(firework2.fw_id)
    # links_dict[firework2.fw_id] = []
    # print(links_dict)

    return Workflow(fireworks, links_dict)

    # return Workflow(fireworks)


def run_lammps_recipe(recipe=LAMMPS_RECIPE_DEFAULT,
                      recipe_settings=RECIPE_SETTING_DEFAULT,
                      working_dir=None,
                      **kwargs):
    """"""
    if not working_dir:
        working_dir=os.getcwd()

    fireworks = []

    links_dict = {}

    for index, step in enumerate(recipe):
        cur_step_dir = os.path.join(working_dir, step[0])
        os.makedirs(cur_step_dir, exist_ok= True)

        if step[1][0] == "template_filename":
            if step[1][1] == "emin_gaff":
                cur_setting = EMIN_SETTING_DEFAULTS.update(recipe_settings[index])
            elif step[1][1] == "npt":
                cur_setting = NPT_SETTING_DEFAULTS.update(recipe_settings[index])
            elif step[1][1] == "nvt":
                cur_setting = NVT_SETTING_DEFAULTS.update(recipe_settings[index])
            else:
                cur_setting = recipe_settings[index]
            cur_firework = RunLammpsFW(working_dir=cur_step_dir,
                                       template_filename = step[1][1],
                                       control_settings = cur_setting,
                                       **kwargs)
        elif step[1][0] == "template_str":
            cur_setting = recipe_settings[index]
            cur_firework = RunLammpsFW(working_dir = cur_step_dir,
                                       template_str = step[1][1],
                                       control_settings = cur_setting,
                                       **kwargs)

        fireworks.append(cur_firework)
        links_dict[cur_firework.fw_id] = []

    firework_ids = list(links_dict.keys())

    for index, firework_id in enumerate(firework_ids[:-1]):
        links_dict[firework_id].append(firework_ids[index + 1])

    return Workflow(fireworks, links_dict)


def fluid_workflow(system_species_data,
                   system_mixture_type,
                   box_data,
                   box_data_type = 'cubic',
                   data_file_name = "complex.data",
                   recipe=LAMMPS_RECIPE_DEFAULT,
                   recipe_settings=RECIPE_SETTING_DEFAULT,
                   recipe_qadapter = QADAPTER_RUN_LAMMPS_SPEC_DEFAULTS,
                   working_dir=None,
                   **kwargs):
    """"""
    if not working_dir:
        working_dir = os.getcwd()

    fireworks = []
    if system_mixture_type == "concentration":
        system_mixture_data = {"Solutes": {},
                               "Solvents": {}}
    elif system_mixture_type == "number of molecules":
        system_mixture_data = {}

    links_dict = {}

    for species, ff_data in system_species_data.items():
        # Generate ff files in separate directories
        # TODO: add logic for ensuring that directory name is legal or at least usable
        cur_species_dir = os.path.join(working_dir, species)
        os.makedirs(cur_species_dir, exist_ok = True)


        # TODO: add logic for not setting file names that do not exist, but code should work as is
        mol2_filename = "{}.mol2".format(species)
        frcmod_filename = "{}.frcmod".format(species)
        prmtop_filename = "{}.prmtop".format(species)
        inpcrd_filename = "{}.inpcrd".format(species)

        cur_firework = GetFFDictFW(ff_data["molecule"],
                                   ff_data["ff_param_data"],
                                   operation_type = ff_data["ff_param_method"],
                                   label = species,
                                   working_dir = cur_species_dir,
                                   output_filename_a = mol2_filename,
                                   input_filename_p = mol2_filename,
                                   output_filename_p = frcmod_filename,
                                   prmtop_filename = prmtop_filename,
                                   tleap_settings = {"mol2_file_path": os.path.join(cur_species_dir, mol2_filename),
                                                     "frcmod_file_path": os.path.join(cur_species_dir, frcmod_filename),
                                                     "prmtop_file_path": os.path.join(cur_species_dir, prmtop_filename),
                                                     "inpcrd_file_path": os.path.join(cur_species_dir, inpcrd_filename)},
                                   **kwargs)

        fireworks.append(cur_firework)

        links_dict[cur_firework.fw_id] = []
        # print(cur_firework.fw_id)
        if system_mixture_type == "concentration":
            system_mixture_data[ff_data["mol_mixture_type"]][species] = ff_data["mixture_data"]
        elif system_mixture_type == "number of molecules":
            system_mixture_data[species] = ff_data["mixture_data"]

    spec = kwargs.pop("spec", {})
    firework2 = Firework([WriteDataFile(working_dir = working_dir,
                                        data_filename = data_file_name,
                                        system_mixture_data_type = system_mixture_type,
                                        system_mixture_data = system_mixture_data,
                                        system_box_data = box_data,
                                        system_box_data_type = box_data_type,
                                        **{i: j for i, j in kwargs.items() if i in
                                           WriteDataFile.required_params +
                                           WriteDataFile.optional_params})],
                         name = "Write_Data_File_FW",
                         parents = fireworks[-1],
                         spec = spec)

    fireworks.append(firework2)
    links_dict[firework2.fw_id] = []

    # print("FIREWORKS:")
    # for firework in fireworks:
    #     # print(firework.fw_id)
    #     # print(firework.fw_id % len(fireworks) + 1)
    #     print(firework.tasks)

    # firework_ids = list(links_dict.keys())
    #
    # for index, firework_id in enumerate(firework_ids[:-1]):
    #     links_dict[firework_id].append(firework_ids[index + 1])
    #     # links_dict[firework_id] += firework_ids[index + 1:]
    #     # links_dict[firework_id].append(firework2.fw_id)
    # links_dict[firework_ids[-1]].append(firework2.fw_id)
    # # links_dict[firework2.fw_id] = []
    # # print(links_dict)

    ### run LAMMPS part
    if not working_dir:
        working_dir=os.getcwd()

    run_fireworks = []

    # links_dict = {}

    for index, step in enumerate(recipe):
        cur_step_dir = os.path.join(working_dir, step[0])
        os.makedirs(cur_step_dir, exist_ok= True)
        print(step[1][0])
        cur_qadapter_spec = QADAPTER_RUN_LAMMPS_SPEC_DEFAULTS[index].copy()
        cur_qadapter_spec.update(recipe_qadapter[index])
        if step[1][0] == "template_filename":
            if step[1][1] == "emin_gaff":
                cur_setting = EMIN_SETTING_DEFAULTS.copy()
            elif step[1][1] == "npt":
                cur_setting = NPT_SETTING_DEFAULTS.copy()
            elif step[1][1] == "nvt":
                cur_setting = NVT_SETTING_DEFAULTS.copy()
            else:
                cur_setting = recipe_settings[index]
            cur_setting.update(recipe_settings[index])
            print('filename')
            cur_firework = RunLammpsFW(working_dir=cur_step_dir,
                                       template_filename = step[1][1],
                                       control_settings = cur_setting,
                                       spec = {"_queueadapter": cur_qadapter_spec},
                                       **kwargs)
        elif step[1][0] == "template_str":
            cur_setting = recipe_settings[index]
            print('string')
            cur_firework = RunLammpsFW(working_dir = cur_step_dir,
                                       template_str = step[1][1],
                                       control_settings = cur_setting,
                                       spec = {"_queueadapter": cur_qadapter_spec},
                                       **kwargs)

        run_fireworks.append(cur_firework)
        links_dict[cur_firework.fw_id] = []

    print("before lammps run fireworks")
    for fw in fireworks:
        print(fw.name, fw.fw_id)
    print("links dict")
    print(links_dict)
    print("\n")

    fireworks += run_fireworks
    firework_ids = list(links_dict.keys())

    print("after lammps run fireworks")
    for fw in fireworks:
        print(fw.name, fw.fw_id)
    print("links dict")
    print(links_dict)
    print("\n")

    for index, firework_id in enumerate(firework_ids[:-1]):
        links_dict[firework_id].append(firework_ids[index + 1])

    print("after reordering of firework ids")
    for fw in fireworks:
        print(fw.name, fw.fw_id)
    print("links dict")
    print(links_dict)

    return Workflow(fireworks, links_dict)

if __name__ == "__main__":
    emin_settings = {"data_file_name": "../complex.data",
                     "restart_final_name": "restart.emin.restart"}
    print(emin_settings)
    print(EMIN_SETTING_DEFAULTS)
    new_settings = EMIN_SETTING_DEFAULTS.copy()
    new_settings.update(emin_settings)
    print(new_settings)
    print(EMIN_SETTING_DEFAULTS)