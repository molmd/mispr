import os

from collections import OrderedDict

from fireworks import Workflow, LaunchPad, FiretaskBase, explicit_serialize
from fireworks.core.rocket_launcher import rapidfire

from pymatgen.io.gaussian import GaussianOutput
from pymatgen.io.ambertools import PrmtopParser
from pymatgen.core.structure import Molecule

from mispr.lammps.workflows.base import lammps_workflow
from mispr.gaussian.utilities.metadata import get_mol_formula


@explicit_serialize
class PrintFW(FiretaskBase):
    """
    Firetask for confirming that modspec works as intended in ProcessPrmtop firetask
    """

    def run_task(self, fw_spec):
        print(str(fw_spec["system_force_field_dict"]))


if __name__ == "__main__":

    # set up the LaunchPad and reset it
    launchpad = LaunchPad(
        host="mongodb://superuser:idlewide@localhost:27017/fireworks?authSource=admin",
        uri_mode=True,
    )
    launchpad.reset("", require_password=False)

    working_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "test_files", "data"
    )

    # Preparing inputs for ambertools for species not using existing ff models
    # Phen_type parameter prep
    Phen_type = "dhps"
    Phen_type_gaussian_output = GaussianOutput(
        os.path.join(working_dir, Phen_type + ".out")
    )
    Phen_type_molecule = Phen_type_gaussian_output.structures[-1]
    Phen_type_molecule.set_charge_and_spin(
        Phen_type_gaussian_output.charge, Phen_type_gaussian_output.spin_multiplicity
    )
    Phen_type_label = get_mol_formula(Phen_type_molecule)
    Phen_type_param_dict = PrmtopParser(
        os.path.join(working_dir, Phen_type_label + ".prmtop"), Phen_type_molecule, ""
    ).to_dict()
    # Oh parameter prep
    Oh_gaussian_output = GaussianOutput(os.path.join(working_dir, "oh.out"))
    Oh_molecule = Oh_gaussian_output.structures[-1]
    Oh_molecule.set_charge_and_spin(
        Oh_gaussian_output.charge, Oh_gaussian_output.spin_multiplicity
    )
    Oh_label = get_mol_formula(Oh_molecule)
    Oh_param_dict = PrmtopParser(
        os.path.join(working_dir, Oh_label + ".prmtop"), Oh_molecule, ""
    ).to_dict()

    # Preparing initial ff dict for species using existing ff models
    # Spce parameter prep
    Spce_molecule = Molecule.from_file(os.path.join(working_dir, "SPC_E.pdb"))
    Spce_label = get_mol_formula(Spce_molecule)
    Spce_param_dict = {
        "Molecule": Spce_molecule,
        "Labels": ["ow", "hw", "hw"],
        "Masses": OrderedDict({"ow": 16.000, "hw": 1.008}),
        "Nonbond": [[0.155394259, 3.16555789], [0.0, 0.0]],
        "Bonds": [{"coeffs": [553.0, 1], "types": [("ow", "hw")]}],
        "Angles": [{"coeffs": [100.0, 109.47], "types": [("hw", "ow", "hw")]}],
        "Dihedrals": [],
        "Impropers": [],
        "Improper Topologies": None,
        "Charges": [-0.8476, 0.4238, 0.4238],
    }
    # Na parameter prep
    Na_molecule = Molecule.from_file(os.path.join(working_dir, "Na.pdb"))
    Na_molecule.set_charge_and_spin(1)
    Na_label = get_mol_formula(Na_molecule)
    Na_param_dict = {
        "Molecule": Na_molecule,
        "Labels": ["na"],
        "Masses": OrderedDict({"na": 22.99}),
        "Nonbond": [[0.02639002, 2.590733472]],  # from frcmod.ions1lm_126_spce (2015)
        #                     'Nonbond': [[0.3526418, 2.159538]], # from frcmod.ionsjc_tip3p (2008)
        "Bonds": [],
        "Angles": [],
        "Dihedrals": [],
        "Impropers": [],
        "Improper Topologies": None,
        "Charges": [1.0],
    }

    sys_ff_dict = {
        Phen_type_label: Phen_type_param_dict,
        Oh_label: Oh_param_dict,
        Na_label: Na_param_dict,
        Spce_label: Spce_param_dict,
    }

    spec = {"system_force_field_dict": sys_ff_dict}

    # setting other inputs for creating the data file
    x = 1.0  # Phen_type concentration

    system_mixture_data_type = "concentration"
    system_mixture_data = {
        "Solutes": {
            Phen_type_label: {
                "Initial Molarity": x,
                "Final Molarity": x,
                "Density": 1.25,
                "Molar Weight": 180.21,
            },
            Oh_label: {
                "Initial Molarity": 3 * x + 1,
                "Final Molarity": 1,
                "Density": 2.13,
                "Molar Weight": 17.007,
            },
            Na_label: {
                "Initial Molarity": 3 * x + 1,
                "Final Molarity": 3 * x + 1,
                "Density": 2.13,
                "Molar Weight": 22.990,
            },
        },
        "Solvents": {Spce_label: {"Density": 0.99705, "Molar Weight": 18.015}},
    }

    # other option for system_mixture_data
    # system_mixture_data_type = 'number of molecules'
    # system_mixture_data = {Phen_type_label: 9,
    #                        Spce_label: 438,
    #                        Oh_label: 9,
    #                        Na_label: 38}

    system_box_side_length = 25.0
    position_seed = 51235
    velocity_seed = 36624
    data_filename = f"complex_from_{system_mixture_data_type}.data"
    sys_mix_type = "concentration"
    # sys_mix_type = "number of molecules"
    sys_species_data = {
        Phen_type_label: {
            "molecule": Phen_type_molecule,
            "ff_param_method": "get_from_prmtop",
            "ff_param_data": os.path.join(working_dir, Phen_type_label + ".prmtop"),
            "save_ff_to_file": True,
            "save_ff_to_db": True,
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
        Spce_label: {
            "molecule": Spce_molecule,
            "ff_param_method": "get_from_dict",
            "ff_param_data": Spce_param_dict,
            "save_ff_to_file": True,
            "mol_mixture_type": "Solvents",
            "mixture_data": {"Density": 0.99705, "Molar Weight": 18.015}
            if sys_mix_type == "concentration"
            else 438,
        },
        Oh_label: {
            "molecule": Oh_molecule,
            "ff_param_method": "get_from_dict",
            "ff_param_data": Oh_param_dict,
            "save_ff_to_file": True,
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
        Na_label: {
            "molecule": Na_molecule,
            "ff_param_method": "get_from_dict",
            "ff_param_data": Na_param_dict,
            "save_ff_to_file": True,
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
        ["nvt_0-20", ["template_filename", "nvt"]],
        ["nvt_20-40", ["template_filename", "nvt"]],
    ]

    v_seed = velocity_seed
    # NEED TO CHANGE FOR EACH NEW PHEN TYPE
    group_def = "\n".join(
        [
            "group phen type 1 2 3 4 5 6 7",
            "group wat type 8 9",
            "group oh type 10 11",
            "group na type 12",
        ]
    )

    shake_group = "wat"
    shake_top = "b 11 a 16"
    timestep = 2
    dump_per = 25000
    compute_def = "\n".join(
        [
            "compute P phen msd com yes",
            "compute W wat msd com yes",
            "compute O oh msd com yes",
            "compute N na msd com yes",
        ]
    )
    ts_comp = "c_P[4] c_W[4] c_O[4] c_N[4]"
    recipe_settings = [
        {
            "data_file_name": "../complex.data",
            "restart_final_name": "restart.emin.restart",
        },
        {
            "restart_file_name": "../emin/restart.emin.restart",
            "group_definitions": group_def,
            "velocity_seed": v_seed,
            "shake_logic": "",
            "shake_group": shake_group,
            "shake_topologies": shake_top,
            "thermo": 1000,
            "timestep": timestep,
            "dump_period": dump_per,
            "restart_final_file_name": "restart.npt.restart",
            "run": 1000000,
        },
        {
            "restart_file_name": "../npt/restart.npt.restart",
            "group_definitions": group_def,
            "shake_logic": "",
            "shake_group": shake_group,
            "shake_topologies": shake_top,
            "temperature_initial": 500.0,
            "temperature_final": 500.0,
            "thermo": 1000,
            "timestep": 1,
            "dump_period": 50000,
            "run": 1000000,
            "restart_final_file_name": "restart.melt_500K.restart",
        },
        {
            "restart_file_name": "../melt/restart.melt_500K.restart",
            "group_definitions": group_def,
            "shake_logic": "",
            "shake_group": shake_group,
            "shake_topologies": shake_top,
            "temperature_initial": 500.0,
            "temperature_final": 298.15,
            "thermo": 1000,
            "timestep": timestep,
            "dump_period": dump_per,
            "run": 1000000,
            "restart_final_file_name": "restart.quench_298-15K.restart",
        },
        {
            "restart_file_name": "../quench/restart.quench_298-15K.restart",
            "group_definitions": group_def,
            "shake_logic": "",
            "shake_group": shake_group,
            "shake_topologies": shake_top,
            "temperature_initial": 298.15,
            "temperature_final": 298.15,
            "thermo": 2,
            "compute_definitions": compute_def,
            "thermo_style_compute": ts_comp,
            "timestep": timestep,
            "dump_period": dump_per,
            "run": 10000000,
            "restart_final_file_name": "restart.nvt_20-ns.restart",
        },
        {
            "restart_file_name": "../nvt_0-20/restart.nvt_20-ns.restart",
            "group_definitions": group_def,
            "shake_logic": "",
            "shake_group": shake_group,
            "shake_topologies": shake_top,
            "temperature_initial": 298.15,
            "temperature_final": 298.15,
            "thermo": 2,
            "compute_definitions": compute_def,
            "thermo_style_compute": ts_comp,
            "timestep": timestep,
            "reset_timestep_logic": "#",
            "dump_period": dump_per,
            "run": 10000000,
            "restart_final_file_name": "restart.nvt_40-ns.restart",
        },
    ]
    nodes = 2
    ntasks_node = 24
    ntasks = nodes * ntasks_node
    mail_user = None
    mail_type = None
    partition = "rajput-24core"

    qadapter = [
        {
            "queue": partition,
            "walltime": "04:00:00",
            "job_name": "emin",
            "nodes": nodes,
            "ntasks": ntasks,
            "ntasks_per_node": ntasks_node,
            "mail-user": mail_user,
            "mail-type": mail_type,
        },
        {
            "queue": partition,
            "walltime": "12:00:00",
            "job_name": "npt",
            "nodes": nodes,
            "ntasks": ntasks,
            "ntasks_per_node": ntasks_node,
            "mail-user": mail_user,
            "mail-type": mail_type,
        },
        {
            "queue": partition,
            "walltime": "06:00:00",
            "job_name": "melt",
            "nodes": nodes,
            "ntasks": ntasks,
            "ntasks_per_node": ntasks_node,
            "mail-user": mail_user,
            "mail-type": mail_type,
        },
        {
            "queue": partition,
            "walltime": "12:00:00",
            "job_name": "quench",
            "nodes": nodes,
            "ntasks": ntasks,
            "ntasks_per_node": ntasks_node,
            "mail-user": mail_user,
            "mail-type": mail_type,
        },
        {
            "queue": partition,
            "walltime": "48:00:00",
            "job_name": "nvt_20",
            "nodes": nodes,
            "ntasks": ntasks,
            "ntasks_per_node": ntasks_node,
            "mail-user": mail_user,
            "mail-type": mail_type,
        },
        {
            "queue": partition,
            "walltime": "48:00:00",
            "job_name": "nvt_40",
            "nodes": nodes,
            "ntasks": ntasks,
            "ntasks_per_node": ntasks_node,
            "mail-user": mail_user,
            "mail-type": mail_type,
        },
    ]

    full_fws, full_links_dict = lammps_workflow(
        system_species_data=sys_species_data,
        system_mixture_type=sys_mix_type,
        box_data=system_box_side_length,
        working_dir=working_dir,
        recipe=recipe,
        recipe_settings=recipe_settings,
        recipe_qadapter=qadapter,
        save_runs_to_file=True,
    )

    print(full_links_dict)
    for firework in full_fws:
        print(firework)

    launchpad.add_wf(Workflow(full_fws, full_links_dict))
    rapidfire(launchpad)
