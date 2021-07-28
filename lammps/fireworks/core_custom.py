# coding: utf-8

# Defines custom fireworks for ambertools and lammps

# import sys
# print(sorted(sys.modules.keys()))

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

def AmbertoolsTasks(mol, **kwargs):

    common_t = [ilftr.RunAntechamber(**{i: j for i, j in kwargs.items() if i in
                                  ilftr.RunAntechamber.required_params + ilftr.RunAntechamber.optional_params}),
                ilftr.RunParmchk(**{i: j for i, j in kwargs.items() if i in
                              ilftr.RunParmchk.required_params + ilftr.RunParmchk.optional_params}),
                ilftwi.WriteTleapScript(**{i: j for i, j in kwargs.items() if i in
                                    ilftwi.WriteTleapScript.required_params + ilftwi.WriteTleapScript.optional_params}),
                ilftr.RunTleap(**{i: j for i, j in kwargs.items() if i in
                            ilftr.RunTleap.required_params + ilftr.RunTleap.optional_params}),
                ilfpo.ProcessPrmtop(molecule = mol,
                              **{i: j for i, j in kwargs.items() if i in
                                 ilfpo.ProcessPrmtop.required_params + ilfpo.ProcessPrmtop.optional_params})]
    # print("common_t:")
    # print(common_t)
    return common_t

class EspToData(fw.Firework):
    def __init__(self,
                 prmtop_filename,
                 file_label,
                 working_dir = None,):
        pass

class GetFFDictFW(fw.Firework):
    def __init__(self,
                 mol,
                 data,
                 operation_type="get_from_esp",
                 db=None,
                 label = "",
                 name="get_ff_dict",
                 parents=None,
                 working_dir=None,
                 **kwargs):
        tasks = []
        working_dir = working_dir or os.getcwd()


        if not label:
            label = iguu.get_mol_formula(mol)

        # print(operation_type)

        if operation_type == "get_from_esp":
            if not isinstance(data, str):
                raise TypeError('data must be a str of the path to the esp file'
                                'for operation_type="get_from_esp"')
            tasks += AmbertoolsTasks(mol,
                                     db=db,
                                     working_dir=working_dir,
                                     input_filename_a=data,
                                     unique_molecule_name=label,
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
                raise TypeError('data must be a str of the path to the prmtop file'
                                'for operation_type="get_from_prmtop"')
            tasks.append(ilfpo.ProcessPrmtop(
                            molecule=mol,
                            working_dir=working_dir,
                            db=db,
                            prmtop_path=data,
                            unique_molecule_name=label,
                            **{i: j for i, j in kwargs.items() if i in
                                        ilfpo.ProcessPrmtop.required_params +
                                        ilfpo.ProcessPrmtop.optional_params}))

        elif operation_type == "get_from_dict":
            if not isinstance(data, dict):
                raise TypeError('data must be a dict for \
                                operation_type="get_from_dict"')
            tasks.append(ilftwi.LabelFFDict(
                            mol=mol,
                            unlabeled_dict=data,
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

class WriteLammpsInputs(fw.Firework):
    def __init__(self,
                 system_mixture_data_type,
                 system_mixture_data,
                 system_box_data,
                 system_box_data_type,
                 control_template,
                 control_settings,
                 working_dir=None,
                 data_filename="complex.data",
                 control_file="control.lammpsin",
                 **kwargs):
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


class RunLammpsFW(fw.Firework):
    def __init__(self,
                 control_file=None,
                 db=None,
                 name="run_lammps",
                 parents=None,
                 working_dir=None,
                 **kwargs):
        """"""
        tasks = []
        working_dir = working_dir or os.getcwd()

        if not control_file:
            tasks.append(ilftwi.WriteControlFile(working_dir=working_dir,
                                                 db=db,
                                          **{i: j for i, j in kwargs.items() if i in
                                             ilftwi.WriteControlFile.required_params + ilftwi.WriteControlFile.optional_params}))

        tasks.append(ilftr.RunLammps(working_dir = working_dir,
                               **{i: j for i, j in kwargs.items() if i in
                                  ilftr.RunLammps.required_params + ilftr.RunLammps.optional_params}))

        super(RunLammpsFW, self).__init__(tasks,
                                          parents=parents,
                                          name=name,
                                          **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS})



class RunAnalysisFW(fw.Firework):
    def __init__(self,
                 property,
                 name="run_analysis",
                 parents=None,
                 working_dir=None,
                 **kwargs):
        """"""
        tasks = []
        working_dir = working_dir or os.getcwd()

        if property == "diffusion":
            msd_kwargs = {i: j for i, j in kwargs.items() if i in ilfpo.GetMSD.required_params + ilfpo.GetMSD.optional_params}
            msd_kwargs.update({"working_dir": working_dir})
            print(msd_kwargs)
            tasks.append(ilfpo.GetMSD(**msd_kwargs))
            diff_kwargs = {i: j for i, j in kwargs.items() if i in ilfpo.CalcDiff.required_params + ilfpo.CalcDiff.optional_params}
            diff_kwargs.update({"working_dir": working_dir})
            tasks.append(ilfpo.CalcDiff(**diff_kwargs))

        elif property == "rdf":
            rdf_kwargs = {i: j for i, j in kwargs.items() if i in inspect.getfullargspec(alsrc.calc_atomic_rdf).args + inspect.getfullargspec(alsrc.calc_molecular_rdf).args + ilfpo.GetRDF.required_params + ilfpo.GetRDF.optional_params}
            rdf_kwargs.update({"working_dir": working_dir})
            tasks.append(ilfpo.GetRDF(**rdf_kwargs))

        super(RunAnalysisFW, self).__init__(tasks,
                                            parents=parents,
                                            name=name,
                                            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS})
#
# fws = []
#
# if file_label:
#     file_label = file_label
# elif molecule:
#     file_label = iguu.get_mol_formula(molecule)
# else not molecule:
#     file_label = esp_file_name.split('.')[0]
#
#
# mol2_filename = f"{file_label}.mol2"
# frcmod_filename = f"{file_label}.frcmod"
# prmtop_filename = f"{file_label}.prmtop"
# inpcrd_filename = f"{file_label}.inpcrd"
#
# antechamber_fw = ilftr.RunAntechamber()

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

    launchpad = LaunchPad(host="mongodb+srv://mbliss01:idlewide@gettingstarted.dt0sv.mongodb.net/fireworks", uri_mode=True)
    launchpad.reset('', require_password=False)

    # firework = GetFFDictFW(dhps_mol,
    #                        esp_file_path,
    #                        working_dir=working_dir,
    #                        output_filename_a=mol2_file,
    #                        prmtop_filename=mol2_file.split('.')[0] + ".prmtop")

    working_dir = "/Users/matt/Documents/GitHub/infrastructure/infrastructure/lammps/tests/test_files/analysis"
    log_file_name = "log.lammps"
    firework = RunAnalysisFW("diffusion",
                             msd_method = "from_log",
                             file_pattern = log_file_name,
                             working_dir = working_dir,
                             kwargs = {"dt": 2})

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
