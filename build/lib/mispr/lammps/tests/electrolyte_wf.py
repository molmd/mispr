import os

from collections import OrderedDict

import numpy as np

from fireworks import LaunchPad
from fireworks.core.fworker import FWorker
from fireworks.queue.queue_launcher import rapidfire
from fireworks.user_objects.queue_adapters.common_adapter import CommonAdapter

from pymatgen.io.gaussian import GaussianOutput
from pymatgen.core.structure import Molecule

from mispr.lammps.workflows.base import lammps_workflow
from mispr.gaussian.utilities.metadata import get_mol_formula

if __name__ == "__main__":

    # set up the LaunchPad and reset it
    launchpad = LaunchPad(
        host="mongodb://superuser:idlewide@10.10.1.100:27017/fireworks?authSource=admin",
        uri_mode=True,
    )
    launchpad.reset("", require_password=False)

    # set up FWorker and QueueAdapterBase
    worker = FWorker.from_file("my_fworker.yaml", f_format="yaml")
    qadapter = CommonAdapter.from_file("my_qadapter.yaml")

    working_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_files")

    # set molecule objects (and labels) for all species
    dhps_gout = GaussianOutput(os.path.join(working_dir, "dhps.out"))
    dhps_mol = dhps_gout.structures[-1]
    dhps_mol.set_charge_and_spin(dhps_gout.charge, dhps_gout.spin_multiplicity)
    dhps_label = get_mol_formula(dhps_mol)

    water_mol = Molecule.from_file(os.path.join(working_dir, "SPC_E.pdb"))
    water_mol.set_charge_and_spin(0, 1)
    water_label = get_mol_formula(water_mol)

    oh_gout = GaussianOutput(os.path.join(working_dir, "oh.out"))
    oh_mol = oh_gout.structures[-1]
    oh_mol.set_charge_and_spin(oh_gout.charge, oh_gout.spin_multiplicity)
    oh_label = get_mol_formula(oh_mol)

    na_mol = Molecule.from_file(os.path.join(working_dir, "Na.pdb"))
    na_mol.set_charge_and_spin(1, 1)
    na_label = get_mol_formula(na_mol)

    # set concentration of phenazine derivative
    x = 1

    # set mixture type ('concentration' or 'number of molecules')
    sys_mix_type = "concentration"
    #    sys_mix_type = 'number of molecules'

    # set the side length of the cube
    box_data = 25.0
    box_data_type = "cubic"

    #    box_data = [[0, 23],
    #                [0, 24],
    #                [0, 25]]
    #    box_data_type = 'rectangular'

    # set the info for how to obtain labeled ff_dict for each species
    # as well as mixture data
    system_species_data = {
        dhps_label: {
            "molecule": dhps_mol,
            "ff_param_method": "get_from_esp",
            "ff_param_data": os.path.join(working_dir, "dhps.esp"),
            "mol_mixture_type": "Solutes",
            "mixture_data": {
                "Initial Molarity": x,
                "Final Molarity": x,
                "Density": 1.25,
                "Molar Weight": 180.21,
            }
            if sys_mix_type == "concentration"
            else 8,
        },
        water_label: {
            "molecule": water_mol,
            "ff_param_method": "get_from_dict",
            "ff_param_data": {
                "Molecule": water_mol,
                "Labels": ["ow", "hw", "hw"],
                "Masses": OrderedDict({"ow": 16.000, "hw": 1.008}),
                "Nonbond": [[0.155394259, 3.16555789], [0.0, 0.0]],
                "Bonds": [{"coeffs": [553.0, 1], "types": [("ow", "hw")]}],
                "Angles": [{"coeffs": [100.0, 109.47], "types": [("hw", "ow", "hw")]}],
                "Dihedrals": [],
                "Impropers": [],
                "Improper Topologies": None,
                "Charges": np.asarray([-0.8476, 0.4238, 0.4238]),
            },
            "mol_mixture_type": "Solvents",
            "mixture_data": {"Density": 0.99705, "Molar Weight": 18.015}
            if sys_mix_type == "concentration"
            else 438,
        },
        oh_label: {
            "molecule": oh_mol,
            "ff_param_method": "get_from_esp",
            "ff_param_data": os.path.join(working_dir, "oh.esp"),
            "mol_mixture_type": "Solutes",
            "mixture_data": {
                "Initial Molarity": 3 * x + 1,
                "Final Molarity": 1,
                "Density": 2.13,
                "Molar Weight": 17.007,
            }
            if sys_mix_type == "concentration"
            else 9,
        },
        na_label: {
            "molecule": na_mol,
            "ff_param_method": "get_from_dict",
            "ff_param_data": {
                "Molecule": na_mol,
                "Labels": ["na"],
                "Masses": OrderedDict({"na": 22.99}),
                "Nonbond": [
                    [0.02639002, 2.590733472]
                ],  # from frcmod.ions1lm_126_spce (2015)
                # 'Nonbond': [[0.3526418, 2.159538]], # from frcmod.ionsjc_tip3p (2008)
                "Bonds": [],
                "Angles": [],
                "Dihedrals": [],
                "Impropers": [],
                "Improper Topologies": None,
                "Charges": np.asarray([1.0]),
            },
            "mol_mixture_type": "Solutes",
            "mixture_data": {
                "Initial Molarity": 3 * x + 1,
                "Final Molarity": 3 * x + 1,
                "Density": 2.13,
                "Molar Weight": 22.990,
            }
            if sys_mix_type == "concentration"
            else 38,
        },
    }

    recipe = [
        ["emin", ["template_filename", "emin_gaff"]],
        ["npt", ["template_filename", "npt"]],
        ["melt", ["template_filename", "nvt"]],
        ["quench", ["template_filename", "nvt"]],
        ["nvt", ["template_filename", "nvt"]],
    ]

    recipe_settings = [
        {
            "data_file_name": "../complex.data",
            "restart_final_name": "restart.emin.restart",
        },
        {
            "restart_file_name": "../emin/restart.emin.restart",
            "restart_final_file_name": "restart.npt.restart",
            "run": 10,
        },
        {
            "restart_file_name": "../npt/restart.npt.restart",
            "temperature_initial": 500.0,
            "temperature_final": 500.0,
            "run": 10,
            "restart_final_file_name": "restart.melt_500K.restart",
        },
        {
            "restart_file_name": "../melt/restart.melt_500K.restart",
            "temperature_initial": 500.0,
            "temperature_final": 298.15,
            "run": 10,
            "restart_final_file_name": "restart.quench_298-15K.restart",
        },
        {
            "restart_file_name": "../quench/restart.quench_298-15K.restart",
            "temperature_initial": 298.15,
            "temperature_final": 298.15,
            "run": 10,
            "restart_final_file_name": "restart.nvt_5-ns.restart",
        },
    ]

    # create workflow for running electrolyte workflow
    workflow = lammps_workflow(
        system_species_data,
        sys_mix_type,
        box_data,
        box_data_type=box_data_type,
        recipe=recipe,
        recipe_settings=recipe_settings,
        working_dir=working_dir,
    )

    # store workflow and launch it using queue
    launchpad.add_wf(workflow)
    rapidfire(launchpad, worker, qadapter, reserve=True)
