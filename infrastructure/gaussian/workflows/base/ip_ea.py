import os
import logging
from queue import Queue

from copy import deepcopy

from fireworks import Firework, Workflow

from infrastructure.gaussian.utils.utils import \
    recursive_relative_to_absolute_path, add_solvent_inputs, \
    check_solvent_inputs
from infrastructure.gaussian.firetasks.parse_outputs import IPEAtoDB
from infrastructure.gaussian.workflows.base.core import common_fw, \
    WORKFLOW_KWARGS

logger = logging.getLogger(__name__)


class Node:
    def __init__(self, state: str, phase: str, num_electrons: int,
                 mol_operation_type=None, mol=None, skip_opt_freq=None,
                 parent: 'Node' = None, ref_charge=None,
                 branch_cation_from_anion: bool = False, h_index: list = None):
        self.phase = phase
        self.state = state
        self.parent = parent

        add_charge = 1 if state.lower() == 'cation' \
            else -1 if state.lower() == 'anion' else 0
        add_charge *= num_electrons
        if parent is None:
            self.added_e = 0
            self.added_h = 0
            self.charge = ref_charge
        else:
            if add_charge < 0 or not branch_cation_from_anion:
                self.added_e = self.parent.added_e - add_charge
                self.added_h = self.parent.added_h
            else:
                self.added_e = self.parent.added_e
                self.added_h = self.parent.added_h + add_charge
            self.charge = self.parent.charge + add_charge

        self.gout_key = \
            f"{self.phase.lower()}_{self.added_e}_{self.added_h}"

        self.fireworks = None
        self.children_nodes = []

        self.mol = None
        if parent is None:
            assert mol_operation_type is not None, \
                "if parent is None, a mol_operation type should be given"
            assert mol is not None, "if parent is None, mol should be given"
            self.mol_operation_type = mol_operation_type
            self.mol = mol
            self.skip = skip_opt_freq
            self.process_mol_func = True
            self.from_fw_spec = False
            self.dir_head = None
        else:
            if branch_cation_from_anion and add_charge > 0:
                h_atom = "[H]"
                self.mol_operation_type = "link_molecules"
                self.mol = {"operation_type":
                                ["get_from_run_dict", "get_from_str"],
                            "mol": [self.parent.gout_key, h_atom],
                            "index": [h_index[self.parent.added_h], 0],
                            "bond_order": 1}
                if num_electrons > 1:
                    for i in range(1, num_electrons):
                        self.mol = {"operation_type":
                                        ["link_molecules", "get_from_str"],
                                    "mol": [self.mol, h_atom],
                                    "index":
                                        [h_index[self.parent.added_h + i], 0],
                                    "bond_order": 1}
            else:
                self.mol_operation_type = 'get_from_run_dict'
                self.mol = self.parent.gout_key
            self.process_mol_func = False
            self.skip = False
            self.from_fw_spec = True
            self.dir_head = self.parent.dir_head
        self.mol_name = None

    def create_fireworks(self, opt_gaussian_inputs, freq_gaussian_inputs,
                         solvent_gaussian_inputs, solvent_properties,
                         working_dir, db, cart_coords,
                         branch_cation_from_anion, **kwargs):
        if "mol_name" in kwargs:
            self.mol_name = kwargs["mol_name"]
            self.dir_head = kwargs["mol_name"]
            kwargs.pop("mol_name")
        if self.mol_name:
            self.mol_name = self.mol_name + f'_{self.phase.lower()}'
        if self.parent is None:
            if "process_mol_func" in kwargs:
                self.process_mol_func = kwargs["process_mol_func"]
                kwargs.pop("process_mol_func")
        else:
            if "process_mol_func" in kwargs:
                kwargs.pop("process_mol_func")
            if not self.mol_name:
                self.mol_name = "{}_{}".format(self.dir_head,
                                               self.phase.lower())

        opt_gins = deepcopy(opt_gaussian_inputs)
        freq_gins = deepcopy(freq_gaussian_inputs)

        gins = [opt_gins, freq_gins]
        check_solvent_inputs(gins)

        if self.phase.lower() == 'solution':
            opt_gins, freq_gins = add_solvent_inputs(gins,
                                                     solvent_gaussian_inputs,
                                                     solvent_properties)
        dir_structure = [self.phase]
        sec_dir_name = f'{self.added_e}e'
        if branch_cation_from_anion:
            sec_dir_name = f'{sec_dir_name}{self.added_h}h'
        dir_structure.append(sec_dir_name)

        local_mol, label, local_fws = \
            common_fw(mol_operation_type=self.mol_operation_type,
                      mol=self.mol,
                      working_dir=working_dir,
                      db=db,
                      opt_gaussian_inputs=opt_gins,
                      freq_gaussian_inputs=freq_gins,
                      cart_coords=cart_coords,
                      oxidation_states=None,
                      process_mol_func=self.process_mol_func,
                      mol_name=self.mol_name,
                      dir_head=self.dir_head,
                      gout_key=self.gout_key,
                      skip_opt_freq=self.skip,
                      dir_structure=dir_structure,
                      from_fw_spec=self.from_fw_spec,
                      charge=self.charge,
                      str_type="smi",
                      **kwargs)
        if self.parent is None and self.dir_head is None:
            self.dir_head = label
        self.fireworks = local_fws
        return

    def branch(self, branching_states, branching_phases, num_of_electrons,
               branch_cation_from_anion, h_index):
        if self.state == 'cation':
            branching_states = [i for i in branching_states if i != 'anion']
        if self.state == 'anion':
            if not branch_cation_from_anion:
                branching_states = [i for i in branching_states if
                                    i != 'cation']
        for state in branching_states:
            self.children_nodes.append(Node(
                state, self.phase, num_of_electrons, parent=self,
                branch_cation_from_anion=branch_cation_from_anion,
                h_index=h_index))
        if self.phase in branching_phases:
            branching_phases.remove(self.phase)
        for phase in branching_phases:
            self.children_nodes.append(Node(
                self.state, phase, 0, parent=self,
                branch_cation_from_anion=branch_cation_from_anion))
        return self.children_nodes


def get_ip_ea(mol_operation_type,
              mol,
              ref_charge,
              single_step=False,
              branch_cation_on_anion=False,
              h_index=None,
              num_electrons=1,
              solvent_gaussian_inputs=None,
              solvent_properties=None,
              states=None,
              phases=None,
              electrode_potentials=None,
              gibbs_elec=-0.001378786,
              db=None,
              name="ip_ea_calculation",
              working_dir=None,
              opt_gaussian_inputs=None,
              freq_gaussian_inputs=None,
              cart_coords=True,
              skip_opt_freq=False,
              **kwargs):
    fws = []
    fireworks_dict = {}
    links_dict = {}
    parents_dict = {}
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)
    mol_object = deepcopy(mol)

    if ref_charge != opt_gaussian_inputs.get("charge", ref_charge):
        raise Exception("The provided reference charge is not consistent with "
                        "the one found in the gaussian input parameters.")

    if states is None:
        states = ["cation", "anion"]
    if phases is None:
        phases = ["gas", "solution"]

    for state in states:
        assert states, 'states list is empty'
        if state.lower() not in ["cation", "anion"]:
            raise ValueError("The provided states are not supported. Supported"
                             " ones are reference, cation, and/or anion.")
    for phase in phases:
        assert phases, 'phases list is empty'
        if phase.lower() not in ["gas", "solution"]:
            raise ValueError("The provided phases are not supported. Supported"
                             " ones are gas and/or solution.")

    if branch_cation_on_anion:
        assert h_index is not None, \
            "index at which to attach hydrogen atom should be provided as input"
        assert len(h_index) == num_electrons, \
            "number of indices at which to attach hydrogen atoms should be " \
            "consistent with number of transfer steps"

    root_node = Node('reference', 'gas' if 'gas' in phases else 'solution',
                     0, mol_operation_type, mol_object,
                     skip_opt_freq=skip_opt_freq, ref_charge=ref_charge,
                     branch_cation_from_anion=branch_cation_on_anion)
    solved_nodes = []
    active_nodes = Queue()
    active_nodes.put(root_node)
    while not active_nodes.empty():
        current_node = active_nodes.get()
        current_node.create_fireworks(
            opt_gaussian_inputs, freq_gaussian_inputs, solvent_gaussian_inputs,
            solvent_properties, working_dir, db, cart_coords,
            branch_cation_on_anion, **kwargs)
        fws += current_node.fireworks
        fireworks_dict[current_node.gout_key] = current_node.fireworks
        if current_node.parent is not None:
            links_dict[current_node.parent.fireworks[-1]] = \
                links_dict.get(current_node.parent.fireworks[-1], []) + \
                [current_node.fireworks[0]]
            parents_dict[current_node.gout_key] = current_node.parent.gout_key

        if single_step:
            addition_electrons = num_electrons
        else:
            addition_electrons = 1
        if abs(current_node.added_e) <= num_electrons or \
                abs(current_node.added_h) <= num_electrons:
            if 'gas' in phases:
                if current_node.phase == 'solution':
                    # no branching is done on solution is gas is in the phases
                    branching_phases = []
                    branching_states = []
                else:
                    branching_phases = ['solution']
                    branching_states = deepcopy(states)
            else:
                branching_phases = []
                branching_states = deepcopy(states)

            if abs(current_node.added_e) == num_electrons:
                if branch_cation_on_anion and \
                        current_node.added_h < num_electrons:
                    branching_states = \
                        [i for i in branching_states if i == 'cation']
                else:
                    branching_states = []
            if abs(current_node.added_h) == num_electrons:
                branching_states = []

            branching_phases = [i for i in branching_phases if i in phases]
            children = current_node.branch(branching_states, branching_phases,
                                           addition_electrons,
                                           branch_cation_on_anion,
                                           h_index)
            for child in children:
                active_nodes.put(child)
        solved_nodes.append(current_node)

    gout_keys = [i.gout_key for i in solved_nodes]
    fw_analysis = Firework(
        IPEAtoDB(num_electrons=num_electrons,
                 steps="single" if single_step else "multi",
                 states=states,
                 phases=phases,
                 solvent_gaussian_inputs=solvent_gaussian_inputs,
                 solvent_properties=solvent_properties,
                 electrode_potentials=electrode_potentials,
                 branch_cation_on_anion=branch_cation_on_anion,
                 keys=gout_keys,
                 root_node_key=root_node.gout_key,
                 gibbs_elec=gibbs_elec,
                 db=db,
                 **{i: j for i, j in kwargs.items()
                    if i in IPEAtoDB.required_params +
                    IPEAtoDB.optional_params}),
        parents=fws[:],
        name="{}-{}".format(root_node.dir_head,
                            "ip_ea_analysis"),
        spec={
            '_launch_dir': os.path.join(working_dir, root_node.dir_head,
                                        "analysis")})
    fws.append(fw_analysis)
    return Workflow(fws,
                    name="{}_{}".format(root_node.dir_head, name),
                    links_dict=links_dict,
                    **{i: j for i, j in kwargs.items()
                       if i in WORKFLOW_KWARGS})
