"""Define the binding energy workflow."""

import os
import logging

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.files import recursive_relative_to_absolute_path
from mispr.gaussian.utilities.inputs import handle_gaussian_inputs
from mispr.gaussian.workflows.base.core import common_fw, WORKFLOW_KWARGS
from mispr.gaussian.firetasks.parse_outputs import BindingEnergytoDB

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.4"


logger = logging.getLogger(__name__)


def get_binding_energies(
    mol_operation_type,
    mol,
    index,
    bond_order=1,
    db=None,
    name="binding_energy_calculation",
    working_dir=None,
    opt_gaussian_inputs=None,
    freq_gaussian_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    cart_coords=True,
    oxidation_states=None,
    skips=None,
    **kwargs
):
    """
    Define a workflow for calculating the binding energy between two molecules.

    * **Fireworks 1 & 2**: Optimize the two molecules in parallel.
    * **Firework 3 & 4**: Run a frequency calculation on each molecule.
    * **Firework 5**: Link the two optimized molecules at a binding site and optimize
      the geometry of the generated complex.
    * **Firework 6**: Run a frequency calculation on the optimized complex.
    * **Firework 7**: Calculate the binding energy of the complex and create BE
      document/json file.

    Args:
        mol_operation_type (list): List of strings of the type of molecule operations.
            See process_mol defined in Defines the binding energy
            ``workflow.mispr/gaussian/utilities/mol.py`` for supported operations; e.g.
            ["get_from_mol", "get_from_file"] to get the first molecule from a
            ``Molecule`` object and the second molecule from a file.
        mol (list): List of the source of the two molecules to be processed. Should
            match the order in ``mol_operation_type``; e.g. if ``mol_operation_type`` is
            ["get_from_mol", "get_from_file"], ``mol`` should be
            [``Molecule``, path to molecule file].
        index (list): List of indices of the two sites in the molecules at which they
            are expected to bind; order should match that in ``mol_operation_type``
            and ``mol``.
        bond_order (int, optional): Bond order to calculate the bond length between the
            two sites. Defaults to 1.
        db (str or dict, optional): Database credentials; could be provided as the path
            to the "db.json" file or in the form of a dictionary; if none is provided,
            attempts to get it from the configuration files.
        name (str, optional): name of the workflow; defaults to
            "binding_energy_calculation".
        working_dir (str, optional): Path of the working directory where any required
            input files can be found and output will be created. Defaults to the current
            working directory.
        opt_gaussian_inputs (dict, optional): Dictionary of Gaussian input parameters
            for the optimization step; e.g.:

            .. code-block:: python

                {
                    "functional": "B3LYP",
                    "basis_set": "6-31G(d)",
                    "route_parameters": {"Opt": None},
                    "link0_parameters": {
                        "%chk": "checkpoint.chk",
                        "%mem": "45GB",
                        "%NProcShared": "24"}
                }

            The above default parameters will be used if not specified.
        freq_gaussian_inputs (dict, optional): Dictionary of Gaussian input parameters
            for the frequency step; default parameters will be used if not specified.
        solvent_gaussian_inputs (str, optional): Gaussian input parameters corresponding
            to the implicit solvent model to be used in the ESP calculations, if any;
            e.g.:

            .. code-block:: python

                "(Solvent=TetraHydroFuran)"

            These parameters should only be specified here and not included in the main
            gaussian_inputs dictionary for each job (i.e. ``opt_gaussian_inputs``,
            ``freq_gaussian_inputs``, etc.). Defaults to None.
        solvent_properties (dict, optional): Additional input parameters to be used in
            the ESP calculations and relevant to the solvent model, if any; for example,
            {"EPS":12}. Defaults to None.
        cart_coords (bool, optional): Uses cartesian coordinates in writing Gaussian
            input files if set to ``True``, otherwise uses z-matrix. Defaults to ``True``.
        oxidation_states (dict, optional): Dictionary of oxidation states that can be
            used in setting the charge and spin multiplicity of the molecule; for
            example: {"Li":1, "O":-2}. Defaults to None.
        skips (list, optional): List of lists of jobs to skip for each molecule; e.g.:
            [["opt", "freq"], ["opt"]]; order should be consistent with that in
            ``mol_operation_type`` and ``mol``. Defaults to None.
        kwargs (keyword arguments): Additional kwargs to be passed to the workflow.

    Returns:
        Workflow
    """
    # TODO: test with different charges and spin multiplicities when
    #  deriving molecules
    # TODO: include an option to use free energy instead of SCF energy
    # mol_operation_type = [], mol = [], index = [], skips = [[], []],
    # mol_name = []
    # order of the indices should be consistent with the order of the mols
    # process_mol_func applies to both input molecules, so if set to False,
    # should give mol_name to each molecule, and if set to True, both will
    # take MolFormula even if mol_name is given to either or both
    fws = []
    labels = []
    working_dir = working_dir or os.getcwd()
    gout_keys = ["mol_1", "mol_2", "mol_linked"]
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    gaussian_inputs = handle_gaussian_inputs(
        {"opt": opt_gaussian_inputs, "freq": freq_gaussian_inputs},
        solvent_gaussian_inputs,
        solvent_properties,
    )
    opt_gaussian_inputs = gaussian_inputs["opt"]
    freq_gaussian_inputs = gaussian_inputs["freq"]

    if skips is None:
        skips = [None, None]
    check_result = []
    for i, j in enumerate(skips):
        if j:
            check_result.append(["final_energy"])
        else:
            check_result.append(None)

    parents = []
    for position, [operation, molecule, key, skip, check, molecule_name] in enumerate(
        zip(
            mol_operation_type,
            mol,
            gout_keys[:2],
            skips,
            check_result,
            kwargs.pop("mol_name", [None, None]),
        )
    ):
        _, label, opt_freq_init_fws = common_fw(
            mol_operation_type=operation,
            mol=molecule,
            working_dir=working_dir,
            db=db,
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            cart_coords=cart_coords,
            oxidation_states=oxidation_states,
            gout_key=key,
            skips=skip,
            mol_name=molecule_name,
            check_result=check,
            **kwargs
        )
        fws += opt_freq_init_fws
        parents.append(len(fws))
        labels.append(label)

    final_mol_label = "{}_{}".format(labels[0], labels[1])

    kwargs.pop("process_mol_func", False)

    _, _, opt_freq_final_fws = common_fw(
        mol_operation_type="link_molecules",
        mol={
            "operation_type": ["get_from_run_dict", "get_from_run_dict"],
            "mol": gout_keys[:2],
            "index": index,
            "bond_order": bond_order,
        },
        working_dir=working_dir,
        db=db,
        filename=final_mol_label,
        opt_gaussian_inputs=opt_gaussian_inputs,
        freq_gaussian_inputs=freq_gaussian_inputs,
        cart_coords=cart_coords,
        oxidation_states=oxidation_states,
        process_mol_func=False,
        mol_name=final_mol_label,
        gout_key=gout_keys[-1],
        from_fw_spec=True,
        **kwargs
    )
    fws += opt_freq_final_fws
    links_dict = {fws[i - 1]: fws[-len(opt_freq_final_fws)] for i in parents}
    fw_analysis = Firework(
        BindingEnergytoDB(
            index=index,
            db=db,
            keys=gout_keys,
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            **{
                i: j
                for i, j in kwargs.items()
                if i
                in BindingEnergytoDB.required_params + BindingEnergytoDB.optional_params
            }
        ),
        parents=fws[:],
        name="{}-{}".format(final_mol_label, "binding_energy_analysis"),
        spec={"_launch_dir": os.path.join(working_dir, "analysis")},
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name="{}_{}".format(final_mol_label, name),
        links_dict=links_dict,
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS}
    )
