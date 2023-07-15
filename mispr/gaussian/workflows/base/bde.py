# coding: utf-8


# Defines the bond dissociation energy workflow.

import os

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.utilities.inputs import handle_gaussian_inputs
from mispr.gaussian.utilities.metadata import get_job_name
from mispr.gaussian.fireworks.break_mol import BreakMolFW
from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.gaussian.firetasks.parse_outputs import BDEtoDB

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.3"


def get_bde(
    mol_operation_type,
    mol,
    ref_charge=0,
    fragment_charges=None,
    bonds=None,
    open_rings=False,
    db=None,
    name="bde_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=False,
    visualize=True,
    **kwargs
):
    """
    Defines a dynamic workflow for calculating the bond dissociation
    energy.

        Firework 1: Optimize the principle molecule.
        Firework 2: Run a frequency analysis.
        Firework 3: Break the bonds in the molecule, find all unique
            fragments, and generate new optimization and frequency
            fireworks for each fragment.
        Firework N: Run a BDE calculation for each fragment and create
            a summary BDE document/json file.

    Note 1: Fireworks 1 and 2 are only present if user does not
    request to skip them.

    Note 2: Charges on the fragments in this workflow are assigned as
        (following the method in the BDE analysis module in pymatgen):
        Neutral molecule: [(0, 0), (1, -1), (-1, 1)]
        Molecule with charge -N: [(-N, 0), (-N+1, -N+2), (-N+2, -N+1), (0, -N)]
        Molecule with charge +N: [(N, 0), (N-1, N-2), (N-2, N-1), (0, N)]

    Note 3: If multiple bonds are being broken but one fragment results
        in an error, the workflow will proceed normally and the final
        document will include all bonds except for the one involving the
        failed fragment.

    Args:
        mol_operation_type (str): the type of molecule operation.
            See process_mol defined in mispr/gaussian/utilities/mol.py
            for supported operations.
        mol (Molecule, GaussianOutput, str, dict): source of the
            molecule to be processed. Should match the mol_operation_type.
        ref_charge (int): charge on the principle molecule
        fragment_charges (list): list of additional charges to consider
            on the fragments besides the default ones. If ref_charge is
            -2, by default all fragments will be calculated with a charge
            of 0, -1, and -2. If the user provides fragment_charges=[-3],
            -3 and 1 will be additionally calculated. If the user
            provides fragment_charges=[-2], this will not cause any
            change since they are already calculated by the workflow;
            defaults to None
        bonds (list): list of tuples of the bonds to break; e.g.
            [(0, 1), (1, 2)] will break the bonds between atoms 0 and 1
            and between atoms 1 and 2; if none is specified, will
            attempt to break all bonds; defaults to None
        open_rings (bool): if True, will open rings encountered during
            fragmentation using OpenBabel's local opt.; defaults to False
        db (str or dict): database credentials; could be provided as
            the path to the db.json file or in the form of a dictionary;
            if none is provided, attempts to get it from the
            configuration files
        name (str): name of the workflow; defaults to "bde_calculation"
        working_dir (str): path of the working directory where any
            required input files can be found and output will be created;
            defaults to the current working directory
        opt_gaussian_inputs (dict): dictionary of Gaussian input
            parameters for the optimization step; e.g.:
            {
                "functional": "B3LYP",
                "basis_set": "6-31G(d)",
                "route_parameters": {"Opt": None},
                "link0_parameters": {
                    "%chk": "checkpoint.chk",
                    "%mem": "45GB",
                    "%NProcShared": "24"}
            }
            the above default parameters will be used if not specified
        freq_gaussian_inputs (dict): dictionary of Gaussian input
            parameters for the frequency step; default parameters will
            be used if not specified
        solvent_gaussian_inputs (str): Gaussian input parameters
            corresponding to the implicit solvent model to be used in
            the ESP calculations, if any; e.g.:
            "(Solvent=TetraHydroFuran)"; these parameters should only
            be specified here and not included in the main
            gaussian_inputs dictionary for each job
            (i.e. opt_gaussian_inputs, freq_gaussian_inputs);
            defaults to None
        solvent_properties (dict): additional input parameters to be
            used in the ESP calculations and relevant to the solvent
            model, if any; e.g., {"EPS":12}; defaults to None
        cart_coords (bool): uses cartesian coordinates in writing
            Gaussian input files if set to True,otherwise uses z-matrix;
            defaults to True
        oxidation_states (dict): dictionary of oxidation states that
            can be used in setting the charge and spin multiplicity of
            the molecule; e.g.: {"Li":1, "O":-2}; defaults to None
        skips (list): list of jobs to skip; e.g.: ["opt", "freq"]; only
            applicable to the principle molecule; defaults to None
        visualize (bool): if True, will generate a summary plot of the
            2D structure of the principle molecule with broken bonds
            highlighted in color, along with a bar plot of the
            corresponding BDEs; requires RDKit to be installed for bond
            highlighting; if RDKit is not found, will throw a
            warning and proceed normally
        **kwargs (keyword arguments): additional kwargs to be passed
            to the workflow

    Returns:
        Workflow
    """
    fws = []
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)
    gout_key = "ref_mol"
    gaussian_inputs = handle_gaussian_inputs(
        {"opt": opt_gaussian_inputs, "freq": freq_gaussian_inputs},
        solvent_gaussian_inputs,
        solvent_properties,
    )
    opt_gaussian_inputs = gaussian_inputs["opt"]
    freq_gaussian_inputs = gaussian_inputs["freq"]

    if skips:
        check_result = ["final_energy", "Enthalpy"]
    else:
        check_result = None

    _, label, opt_freq_fws = common_fw(
        mol_operation_type=mol_operation_type,
        mol=mol,
        charge=ref_charge,
        working_dir=working_dir,
        dir_structure=["principle_mol"],
        db=db,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        skips=skips,
        check_result=check_result,
        gout_key=gout_key,
        **kwargs
    )
    fws += opt_freq_fws

    break_fw = BreakMolFW(
        mol=gout_key,
        mol_operation_type="get_from_run_dict",
        from_fw_spec=True,
        bonds=bonds,
        open_rings=open_rings,
        ref_charge=ref_charge,
        fragment_charges=fragment_charges,
        db=db,
        calc_frags=True,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        cart_coords=cart_coords,
        name=get_job_name(label, "breaking"),
        parents=fws[:],
        working_dir=os.path.join(working_dir, label, "fragments"),
        **kwargs
    )
    fws.append(break_fw)

    fw_analysis = Firework(
        BDEtoDB(
            principle_mol_key=gout_key,
            db=db,
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            visualize=visualize,
            **{
                i: j
                for i, j in kwargs.items()
                if i in BDEtoDB.required_params + BDEtoDB.optional_params
            }
        ),
        parents=fws[:],
        name="{}-{}".format(label, "bde_analysis"),
        spec={
            "_launch_dir": os.path.join(working_dir, label, "analysis"),
            "_allow_fizzled_parents": True,
        },
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name="{}_{}".format(label, name),
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS}
    )
