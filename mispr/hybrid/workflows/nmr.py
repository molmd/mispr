# coding: utf-8


# Defines the DFT-MD workflow that extracts solvation structures and computes their
# nuclear magnetic resonances.

import os

from copy import deepcopy

from fireworks import Firework, Workflow

from mispr.hybrid.defaults import (
    NMR_GAUSSIAN_INPUTS,
    OPT_GAUSSIAN_INPUTS,
    FREQ_GAUSSIAN_INPUTS,
)
from mispr.hybrid.workflows.core import run_hybrid_calcs
from mispr.hybrid.firetasks.nmr_from_md import NMRFromMD

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Oct 2021"
__version__ = "0.0.1"


def get_solvation_structures_nmr(
    mol_operation_type,
    mol,
    mol_type,
    mol_data,
    box_data,
    ff_method=None,
    ff_params=None,
    mixture_type="number of molecules",
    db=None,
    name="hybrid_calculation",
    working_dir=None,
    opt_esp_gaussian_inputs=None,
    freq_esp_gaussian_inputs=None,
    esp_gaussian_inputs=None,
    opt_nmr_gaussian_inputs=None,
    freq_nmr_gaussian_inputs=None,
    nmr_gaussian_inputs=None,
    esp_solvent_gaussian_inputs=None,
    esp_solvent_properties=None,
    nmr_solvent_gaussian_inputs=None,
    nmr_solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    box_data_type="cubic",
    data_file_name="data.mixture",
    analysis_list=None,
    analysis_settings=None,
    **kwargs,
):
    """
    Returns a workflow to run hybrid DFT/MD calculations, extract categorize atomic
    clusters around a particle of interest, and run the NMR calculations for the
    clusters that have the highest frequency of occurrence in solution.

    Args:
        mol_operation_type (list of str): operation to perform for each
            molecule composing the liquid solution in order to process the input
            structure format. Length should correspond to the number of molecules/species
            composing the liquid solution. Supported commands:
            1. "get_from_mol" if the input is a pymatgen Molecule object,
            2. "get_from_file" if the input is any file format supported by Openabel and pymatgen,
            3. "get_from_gout_file" if the input is a Gaussian output file,
            4. "get_from_str" if the input is a string,
            5. "get_from_mol_db" if the input is an InChI representation of the molecule to be used to query the database,
            6. "get_from_gout" if the input is a pymatgen GaussianOutput object,
            7. "get_from_run_dict" if the input is a GaussianOutput dictionary,
            8. "get_from_run_id" if the input is a MongoDB document ID to be used to query the database,
            9. "get_from_run_query" if the input is a dictionary with criteria to search the database: e.g. "{'inchi': inchi, 'type': type, 'functional': func, ...}
            10. "derive_molecule": used for deriving a molecule by attaching a functional group at a site and the corresponding mol input should be a dictionary,
                e.g. {'operation_type': <mol_operation_type for the base structure>, 'mol': <base_mol>,'func_grp': func_group_name, ....}
            11."link_molecules": used for linking two structures by forming a bond at specific sites and the corresponding mol input should be a dictionary,
                e.g. {'operation_type': ['get_from_file', 'get_from_mol_db'], 'mol': ['mol1.xyz', 'mol_inchi'], 'index': [3, 5], 'bond_order': 1}
        mol (list): sources of structures making up the liquid solution,
            e.g. file path if mol_operation_type is specified as "get_from_file",
            InChI string if mol_operation_type is specified as "get_from_mol_db", etc.
            Order should match that in mol_operation_type
        mol_type (list of str): type of each structure composing the liquid solution.
            Supported types: "Solvents", "Solutes". Used for calculating the number of
            molecules of each type if this information is not provided.
        mol_data (list of int or list of dict): format depends on mixture_type input.
            If mixture_type is "number of molecules", mol_data should be a list of
            the number of molecules of each type. If mixture_type is "concentration",
            mol_data should be a list of dict, where each dictionary should follow the format:
                {'Initial Molarity': molarity_i,
                 'Final Molarity': molarity_f,
                 'Density': density,
                 'Molar Weight': molar_weight
                }
        box_data (float, int, list (3,2), array (3,2), or LammpsBox): definitions for
            box size. See box_data_type for info how to define this parameter.
        ff_method (list of str): operation to perform for each molecule composing the liquid
            solution in order to process the force field parameters. can be "get_from_esp",
            "get_from_prmtop", or "get_from_dict". Defaults to "get_from_esp" for all molecules.
        ff_params (list of str or dict): sources of the force field parameters for each
            molecule type; type depends on what is specified in the ff_method input;
            if "get_from_esp" is used, the corresponding ff_param should be an
            empty dictionary since the path to the ESP file is automatically detected
            from the ESP calculations that is performed at the beginning of the workflow;
            if "get_from_prmtop" is used, the corresponding ff_param should be the
            path to the prmtop file; if "get_from_dict" is used, the corresponding
            ff_param should be a dictionary, for example:
            {
                "Labels": ["mg"],
                "Masses": OrderedDict({"mg": 24.305}),
                "Nonbond": [[0.8947000005260684, 1.412252647723565]],
                "Bonds": [],
                "Angles": [],
                "Dihedrals": [],
                "Impropers": [],
                "Improper Topologies": None,
                "Charges": np.asarray([2.0]),
            }
        mixture_type (str): "concentration" or "number of molecules"
        db (str or dict): database credentials; could be provided as the path to the
            db.json file or in the form of a dictionary; if none is provided, attempts
            to get it from the configuration files
        name (str): name of the workflow; defaults to "hybrid_calculation"
        working_dir (str): path of the working directory where any required input files
            can be found and output will be created
        opt_esp_gaussian_inputs (dict): dictionary of Gaussian input parameters for the
            optimization step of the ESP workflow; for example:
            {
            "functional": "HF",
            "basis_set": "6-31G(d)",
            "route_parameters": {"Opt": None},
            "link0_parameters": {
                "%chk": "checkpoint.chk",
                "%mem": "45GB",
                "%NProcShared": "28",
            }
            default parameters will be used if not specified
        freq_esp_gaussian_inputs (dict): dictionary of Gaussian input parameters for the
            frequency step of the ESP workflow; default parameters will be used if not specified
        esp_gaussian_inputs (dict): dictionary of Gaussian input parameters for the
            ESP step of the ESP workflow; default parameters will be used if not specified
        opt_nmr_gaussian_inputs (dict): dictionary of Gaussian input parameters for the
            optimization step of the NMR workflow; default parameters will be used if not specified
        freq_nmr_gaussian_inputs (dict): dictionary of Gaussian input parameters for the
            frequency step of the NMR workflow; default parameters will be used if not specified
        nmr_gaussian_inputs (dict): dictionary of Gaussian input parameters for the
            NMR step of the NMR workflow; default parameters will be used if not specified
        esp_solvent_gaussian_inputs (str): Gaussian input parameters corresponding to the
            implicit solvent model to be used in the ESP calculations, if any; for example: "(Solvent=TetraHydroFuran)";
            these parameters should only be specified here and not included in the main
            gaussian_inputs dictionary for each job (i.e. opt_gaussian_inputs, freq_esp_gaussian_inputs, etc.)
        esp_solvent_properties (dict): additional input parameters to be used in the ESP calculations and
            relevant to the solvent model, if any; for example, {"EPS":12}
        nmr_solvent_gaussian_inputs (str): Gaussian input parameters corresponding to the
            implicit solvent model to be used in the NMR calculations, if any; for example: "(Solvent=TetraHydroFuran)";
            these parameters should only be specified here and not included in the main
            gaussian_inputs dictionary for each job (i.e. opt_nmr_gaussian_inputs, freq_nmr_gaussian_inputs, etc.)
        nmr_solvent_properties (dict): additional input parameters to be used in the NMR calculations and
            relevant to the solvent model, if any; for example, {"EPS":12}
        cart_coords (bool): uses cartesian coordinates in writing Gaussian input files if set to True,
            otherwise uses z-matrix; defaults to True
        oxidation_states (dict): dictionary of oxidation states that can be used in setting the charge and
            spin multiplicity of the clusters extracted from MD simulations to be used in the NMR workflow
        skips (list of lists): type of DFT calculation to skip in the ESP workflow for each molecule;
            for example ["opt", "freq"], ["opt"], ["freq"], or [];
        box_data_type (str): Can be one of the following: 'cubic','rectangular', or 'LammpsBox'.
            If 'cubic', box_data must be a float or int; if 'rectangular', box_data must be an array-like
            with size (3,2); if 'LammpsBox', box_data must be a pymatgen.io.lammps.data.LammpsBox object.
            Defaults to 'cubic'.
        data_file_name (str): name of the LAMMPS data file to create and use; defaults to "data.mixture"
        analysis_list (list of str): type of MD analysis to perform after the MD simulations are finished;
            for example: ["diffusion", "rdf", "cn", "clusters"] if user wants to perform diffusion,
            RDF, coordination number, and cluster analysis
        analysis_settings (list of dict): settings of the MD analysis steps; please
            refer to the mdproptools documentation for details of inputs used in the analysis functions;
            order of settings should correspond to the order used in "analysis_list"
        **kwargs keyword arguments): additional kwargs to be passed to the workflow;
            for example: lammps "recipe" and "recipe_settings"; the defaults for these
            inputs are specified in the mispr/lammps/defaults.py

    Returns:
        Workflow
    """
    if not working_dir:
        working_dir = os.getcwd()
    if not analysis_list:
        analysis_list = []
    # analysis_list += ["diffusion", "rdf", "cn", "clusters"]
    # analysis_list = list(set(analysis_list))

    wf = run_hybrid_calcs(
        mol_operation_type,
        mol,
        mol_type,
        mol_data,
        box_data,
        ff_method,
        ff_params,
        mixture_type,
        db,
        name,
        working_dir,
        opt_esp_gaussian_inputs,
        freq_esp_gaussian_inputs,
        esp_gaussian_inputs,
        esp_solvent_gaussian_inputs,
        esp_solvent_properties,
        cart_coords,
        oxidation_states,
        skips,
        box_data_type,
        data_file_name,
        analysis_list,
        analysis_settings,
        **kwargs,
    )

    dft_md_fw_ids = list(wf.id_fw.keys())

    nmr_dir = f"{working_dir}/nmr"
    nmr_fw = Firework(
        NMRFromMD(
            db=db,
            working_dir=nmr_dir,
            opt_gaussian_inputs=opt_nmr_gaussian_inputs or OPT_GAUSSIAN_INPUTS,
            freq_gaussian_inputs=freq_nmr_gaussian_inputs or FREQ_GAUSSIAN_INPUTS,
            nmr_gaussian_inputs=nmr_gaussian_inputs or NMR_GAUSSIAN_INPUTS,
            solvent_gaussian_inputs=nmr_solvent_gaussian_inputs,
            solvent_properties=nmr_solvent_properties,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            additional_kwargs=kwargs,
        ),
        name="nmr_calculation",
        spec={"_launch_dir": nmr_dir},
    )
    wf.append_wf(Workflow.from_Firework(nmr_fw), dft_md_fw_ids[:-1])
    return wf
