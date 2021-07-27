import os

from fireworks.core.firework import Firework, Workflow

from mispr.lammps.fireworks.custom_fw import GetFFDictFW
from mispr.lammps.firetasks.write_inputs import WriteDataFile

__author__ = "Matthew Bliss"
__maintainer__ = "Matthew Bliss"
__email__ = "matthew.bliss@tufts.edu"
__status__ = "Development"
__date__ = "Apr 14, 2020"
__version__ = "0.0.1"


def write_lammps_data(
    system_species_data,
    system_mixture_type,
    box_data,
    box_data_type="cubic",
    data_file_name="complex.data",
    working_dir=None,
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

    fireworks = []
    if system_mixture_type == "concentration":
        system_mixture_data = {"Solutes": {}, "Solvents": {}}
    elif system_mixture_type == "number of molecules":
        system_mixture_data = {}

    links_dict = {}

    for species, ff_data in system_species_data.items():
        mol2_filename = "{}.mol2".format(species)
        frcmod_filename = "{}.frcmod".format(species)
        prmtop_filename = "{}.prmtop".format(species)
        inpcrd_filename = "{}.inpcrd".format(species)

        cur_firework = GetFFDictFW(
            ff_data["molecule"],
            ff_data["ff_param_data"],
            operation_type=ff_data["ff_param_method"],
            label=species,
            working_dir=working_dir,
            output_filename_a=mol2_filename,
            input_filename_p=mol2_filename,
            output_filename_p=frcmod_filename,
            prmtop_filename=prmtop_filename,
            tleap_settings={
                "mol2_file_path": os.path.join(working_dir, mol2_filename),
                "frcmod_file_path": os.path.join(working_dir, frcmod_filename),
                "prmtop_file_path": os.path.join(working_dir, prmtop_filename),
                "inpcrd_file_path": os.path.join(working_dir, inpcrd_filename),
            },
            **kwargs
        )

        fireworks.append(cur_firework)

        links_dict[cur_firework.fw_id] = []
        # print(cur_firework.fw_id)
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
