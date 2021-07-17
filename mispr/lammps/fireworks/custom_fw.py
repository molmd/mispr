# coding: utf-8

# Defines custom fireworks for ambertools and lammps

import os
import pathlib
import logging

from fireworks import Firework, Workflow
from pymatgen.core.structure import Molecule
from infrastructure.lammps.firetasks.write_inputs import WriteDataFile, WriteControlFile, WriteTleapScript,\
    TLEAP_SETTING_DEFAULTS, LabelFFDict
from infrastructure.lammps.firetasks.run import RunLammps, RunAntechamber, RunParmchk, RunTleap
from infrastructure.lammps.firetasks.parse_outputs import ProcessPrmtop
from infrastructure.lammps.utils.utils import add_ff_labels_to_dict
from infrastructure.gaussian.utils.utils import get_mol_formula

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames

def AmbertoolsTasks(mol, **kwargs):

    common_t = [RunAntechamber(**{i: j for i, j in kwargs.items() if i in
                                  RunAntechamber.required_params + RunAntechamber.optional_params}),
                RunParmchk(**{i: j for i, j in kwargs.items() if i in
                              RunParmchk.required_params + RunParmchk.optional_params}),
                WriteTleapScript(**{i: j for i, j in kwargs.items() if i in
                                    WriteTleapScript.required_params + WriteTleapScript.optional_params}),
                RunTleap(**{i: j for i, j in kwargs.items() if i in
                            RunTleap.required_params + RunTleap.optional_params}),
                ProcessPrmtop(molecule = mol,
                              **{i: j for i, j in kwargs.items() if i in
                                 ProcessPrmtop.required_params + ProcessPrmtop.optional_params})]
    # print("common_t:")
    # print(common_t)
    return common_t

class EspToData(Firework):
    def __init__(self,
                 prmtop_filename,
                 file_label,
                 working_dir = None,):
        pass

class GetFFDictFW(Firework):
    def __init__(self,
                 mol,
                 data,
                 operation_type="get_from_esp",
                 label = "",
                 name="get_ff_dict",
                 parents=None,
                 working_dir=None,
                 **kwargs):
        tasks = []
        working_dir = working_dir or os.getcwd()
        if not label:
            label = get_mol_formula(mol)

        # print(operation_type)

        if operation_type == "get_from_esp":
            if not isinstance(data, str):
                raise TypeError('data must be a str of the path to the esp file'
                                'for operation_type="get_from_esp"')
            tasks += AmbertoolsTasks(mol,
                                     working_dir=working_dir,
                                     input_filename_a=data,
                                     **{i: j for i, j in kwargs.items() if i in
                                        RunAntechamber.required_params +
                                        RunAntechamber.optional_params +
                                        RunParmchk.required_params +
                                        RunParmchk.optional_params +
                                        WriteTleapScript.required_params +
                                        WriteTleapScript.optional_params +
                                        RunTleap.required_params +
                                        RunTleap.optional_params +
                                        ProcessPrmtop.required_params +
                                        ProcessPrmtop.optional_params})

        elif operation_type == "get_from_prmtop":
            if not isinstance(data, str):
                raise TypeError('data must be a str of the path to the prmtop file'
                                'for operation_type="get_from_prmtop"')
            tasks.append(ProcessPrmtop(molecule=mol,
                                       working_dir=working_dir,
                                       prmtop_path=data,
                                       **{i: j for i, j in kwargs.items() if i in
                                          ProcessPrmtop.required_params +
                                          ProcessPrmtop.optional_params}))

        elif operation_type == "get_from_dict":
            if not isinstance(data, dict):
                raise TypeError('data must be a dict for operation_type="get_from_dict"')
            tasks.append(LabelFFDict(mol=mol,
                                     unlabeled_dict=data,
                                     working_dir=working_dir,
                                     label=label,
                                     **{i: j for i, j in kwargs.items() if i in
                                        LabelFFDict.required_params +
                                        LabelFFDict.optional_params}))

        # print(tasks)

        super(GetFFDictFW, self).__init__(tasks,
                                          parents=parents,
                                          name=name,
                                          **{i: j for i, j in kwargs.items()
                                             if i in FIREWORK_KWARGS})

class WriteLammpsInputs(Firework):
    def __init__(self,
                 system_mol_data,
                 mixture_type,
                 box_side_length,
                 working_dir=None,
                 data_file="complex.data",
                 control_file="control.lammpsin",
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


#
# fws = []
#
# if file_label:
#     file_label = file_label
# elif molecule:
#     file_label = get_mol_formula(molecule)
# else not molecule:
#     file_label = esp_file_name.split('.')[0]
#
#
# mol2_filename = f"{file_label}.mol2"
# frcmod_filename = f"{file_label}.frcmod"
# prmtop_filename = f"{file_label}.prmtop"
# inpcrd_filename = f"{file_label}.inpcrd"
#
# antechamber_fw = RunAntechamber()

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
    dhps_mol.set_charge_and_spin(dhps_gout.charge,
                                 dhps_gout.spin_multiplicity)

    launchpad = LaunchPad(host="mongodb://superuser:idlewide@localhost:27017/fireworks?authSource=admin", uri_mode=True)
    launchpad.reset('', require_password=False)

    firework = GetFFDictFW(dhps_mol,
                           esp_file_path,
                           working_dir=working_dir,
                           output_filename_a=mol2_file,
                           prmtop_filename=mol2_file.split('.')[0] + ".prmtop")



    launchpad.add_wf(firework)
    launch_rocket(launchpad, FWorker())
