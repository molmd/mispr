import os
import logging

from fireworks import Firework, Workflow

from infrastructure.gaussian.utils.utils import get_mol_formula, \
    recursive_relative_to_absolute_path
from infrastructure.gaussian.firetasks.parse_outputs import BindingEnergytoDB
from infrastructure.gaussian.workflows.base.core import common_fw, \
    WORKFLOW_KWARGS

logger = logging.getLogger(__name__)


def get_binding_energies(mol_operation_type,
                         mol,
                         index,
                         bond_order=1,
                         db=None,
                         name="binding_energy_calculation",
                         working_dir=None,
                         opt_gaussian_inputs=None,
                         freq_gaussian_inputs=None,
                         cart_coords=True,
                         oxidation_states=None,
                         skip_opt_freq=None,
                         **kwargs):
    # TODO: test with different charges and spin multiplicities when
    #  deriving molecules
    # mol_operation_type = [], mol = [], index = []
    # order of the indices should be consistent with the order of the mols
    fws = []
    molecules = []
    working_dir = working_dir or os.getcwd()
    keys = ["mol_1", "mol_2", "mol_linked"]
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    if skip_opt_freq is None:
        skip_opt_freq = [False, False]
    parents = []
    for position, [operation, molecule, key, skip] in \
            enumerate(
                zip(mol_operation_type, mol, keys[:2], skip_opt_freq)):
        mol_object, opt_freq_init_fws = \
            common_fw(mol_operation_type=operation,
                      mol=molecule,
                      working_dir=working_dir,
                      db=db,
                      opt_gaussian_inputs=opt_gaussian_inputs,
                      freq_gaussian_inputs=freq_gaussian_inputs,
                      cart_coords=cart_coords,
                      oxidation_states=oxidation_states,
                      gout_key=key,
                      skip_opt_freq=skip,
                      **kwargs)
        fws += opt_freq_init_fws
        parents.append(len(fws))
        molecules.append(mol_object)

    final_mol_formula = "{}_{}".format(get_mol_formula(molecules[0]),
                                       get_mol_formula(molecules[1]))

    _, opt_freq_final_fws = common_fw(mol_operation_type="link_molecules",
                                      mol={"operation_type": [
                                          "get_from_run_dict",
                                          "get_from_run_dict"],
                                          "mol": keys[:2],
                                          "index": index,
                                          "bond_order": bond_order},
                                      working_dir=working_dir,
                                      db=db,
                                      filename=final_mol_formula,
                                      opt_gaussian_inputs=opt_gaussian_inputs,
                                      freq_gaussian_inputs=freq_gaussian_inputs,
                                      cart_coords=cart_coords,
                                      oxidation_states=oxidation_states,
                                      process_mol_func=False,
                                      mol_name=final_mol_formula,
                                      gout_key=keys[-1],
                                      from_fw_spec=True,
                                      **kwargs)
    fws += opt_freq_final_fws
    links_dict = {fws[i - 1]: fws[-len(opt_freq_final_fws)] for i in parents}
    fw_analysis = Firework(
        BindingEnergytoDB(index=index,
                          db=db,
                          **{i: j for i, j in kwargs.items()
                             if i in BindingEnergytoDB.required_params +
                             BindingEnergytoDB.optional_params}),
        parents=fws[:],
        name="{}-{}".format(final_mol_formula,
                            "binding_energy_analysis"),
        spec={'_launch_dir': os.path.join(working_dir, final_mol_formula,
                                          'Analysis')})
    fws.append(fw_analysis)

    return Workflow(fws,
                    name="{}_{}".format(final_mol_formula, name),
                    links_dict=links_dict,
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})
