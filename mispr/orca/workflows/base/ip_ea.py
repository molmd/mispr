"""Define the redox potentials workflow, using ORCA instead of Gaussian.

Mirrors ``mispr.gaussian.workflows.base.ip_ea.get_ip_ea``'s ``Node``-based tree:
each node is one molecule state (a given charge, in a given phase) and the tree
branches into cation/anion states and gas/solution phases as it grows, so the
same code builds whatever combination of vertical/adiabatic, single/multi-step,
and PCET calculations the caller asks for. ``IPEAtoDB`` (the final analysis
Firetask that converts Gibbs free energies into redox potentials) is reused
unmodified from the Gaussian firetasks -- like ``BDEtoDB``/``BindingEnergytoDB``
elsewhere in this package, it only reads the shared gout_dict schema that
``RunOrca`` also produces, with no Gaussian-specific logic.

Two differences from the Gaussian tree, both driven by how ``RunOrca`` works:

* Implicit solvent is a plain keyword (``CPCM(<solvent>)``) added to whichever
  job is already running, not a separate Gaussian-style PCM route section --
  so unlike the Gaussian workflow, no ``handle_gaussian_inputs`` merging step
  is needed; the solvent dict is just forwarded to ``common_fw``/``OrcaFW``
  alongside the functional/basis inputs.
* "Vertical" (``skips=["opt"]``) still runs a Freq calculation, not a bare
  single point, on the frozen parent geometry -- because IPEAtoDB needs a
  Gibbs free energy for every state, and only a Freq job produces one. This
  matches ``mispr.gaussian.workflows.base.core.common_fw``'s own skips=["opt"]
  branch, which likewise keeps Freq and only drops Opt.
"""

import os
import logging

from copy import deepcopy
from queue import Queue

from fireworks import Firework, Workflow

from mispr.gaussian.utilities.files import (
    bibtex_parser,
    recursive_relative_to_absolute_path,
)
from mispr.gaussian.utilities.mol import process_mol
from mispr.gaussian.firetasks.geo_transformation import ProcessMoleculeInput
from mispr.gaussian.firetasks.parse_outputs import IPEAtoDB
from mispr.orca.firetasks.run_calc import RunOrca
from mispr.orca.fireworks.core import OrcaFW
from mispr.orca.workflows.base.core import common_fw, WORKFLOW_KWARGS

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_8_9"
__version__ = "0.0.5"

logger = logging.getLogger(__name__)


def _to_orca_solvent(solvent_gaussian_inputs):
    """Translate the Gaussian-style solvent string (e.g. "(Solvent=Water)") into the dict RunOrca expects (e.g. {"solvent": "water"}); returns None for gas phase.

    Same translation used in ``mispr.orca.workflows.base.binding_energy``.
    """
    if not solvent_gaussian_inputs:
        return None
    solvent_inputs = [
        i.lower() for i in solvent_gaussian_inputs.strip("()").split(",")
    ]
    solvent_name = next(
        (s.split("=")[1] for s in solvent_inputs if "solvent" in s), "water"
    )
    return {"solvent": solvent_name}


class Node:
    """Generate the Fireworks corresponding to different molecule states in the IP/EA workflow.

    Each molecule state corresponds to a node in the tree. The node is a leaf
    if it is the last node in the tree, otherwise it is a branch. Not meant to
    be instantiated directly. Mirrors
    ``mispr.gaussian.workflows.base.ip_ea.Node``.
    """

    def __init__(
        self,
        state: str,
        phase: str,
        num_electrons: int,
        mol=None,
        skips=None,
        parent: "Node" = None,
        ref_charge=None,
        branch_cation_from_anion: bool = False,
        h_index: list = None,
    ):
        """Set up this node's charge/electron/hydrogen bookkeeping relative to its parent.

        Args:
            state (str): Current state of the molecule: cation or anion.
            phase (str): Current phase of the molecule: gas or solution.
            num_electrons (int): Number of electrons to transfer.
            mol (Molecule, optional): Already-resolved pymatgen Molecule for
                the root node; required if ``parent`` is None, ignored
                otherwise (child nodes chain onto their parent's ORCA run
                instead).
            skips (list, optional): List of jobs to skip; only ever ``None``
                or ``["opt"]`` in this workflow (the latter for vertical
                calculations). Defaults to None.
            parent (Node, optional): Parent node of the current node; None if
                the node corresponds to the initial molecule state.
            ref_charge (int, optional): The initial charge on the molecule;
                only relevant to the parent node.
            branch_cation_from_anion (bool, optional): Whether to add a
                hydrogen atom at the current node; relevant for PCET
                calculations.
            h_index (list, optional): Site indices in the molecule at which to
                attach the hydrogen atoms in the PCET calculations.

        """
        self.phase = phase
        self.state = state
        self.parent = parent

        add_charge = (
            1 if state.lower() == "cation" else -1 if state.lower() == "anion" else 0
        )
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

        self.gout_key = f"{self.phase.lower()}_{self.added_e}_{self.added_h}"

        self.fireworks = None
        self.children_nodes = []
        self.dir_head = None
        self.mol_name = None

        if parent is None:
            if mol is None:
                raise ValueError("if parent is None, mol should be given")
            self.mol = mol
            self.link_mol = None
        else:
            self.dir_head = self.parent.dir_head
            self.mol = None
            if branch_cation_from_anion and add_charge > 0:
                h_atom = "[H]"
                self.link_mol = {
                    "operation_type": ["get_from_run_dict", "get_from_str"],
                    "mol": [self.parent.gout_key, h_atom],
                    "index": [h_index[self.parent.added_h], 0],
                    "bond_order": 1,
                }
                if num_electrons > 1:
                    for i in range(1, num_electrons):
                        self.link_mol = {
                            "operation_type": ["link_molecules", "get_from_str"],
                            "mol": [self.link_mol, h_atom],
                            "index": [h_index[self.parent.added_h + i], 0],
                            "bond_order": 1,
                        }
            else:
                self.link_mol = None
        self.skip = skips

    def create_fireworks(
        self,
        opt_orca_inputs,
        freq_orca_inputs,
        solvent,
        working_dir,
        db,
        branch_cation_from_anion,
        tag,
        orca_settings,
        **kwargs,
    ):
        """Generate the optimization and/or frequency Fireworks corresponding to the current node."""
        if "mol_name" in kwargs:
            self.mol_name = kwargs.pop("mol_name")
            self.dir_head = self.mol_name
        if self.mol_name:
            self.mol_name = f"{self.mol_name}_{self.phase.lower()}"
        elif self.parent is not None:
            self.mol_name = f"{self.dir_head}_{self.phase.lower()}"

        node_solvent = solvent if self.phase.lower() == "solution" else None

        dir_structure = [self.phase]
        sec_dir_name = f"{self.added_e}e"
        if branch_cation_from_anion:
            sec_dir_name = f"{sec_dir_name}{self.added_h}h"
        dir_structure.append(sec_dir_name)

        if self.link_mol is not None:
            # PCET step: attach an H atom to the parent's optimized geometry,
            # then optimize (or, for vertical, just run Freq on the linked-but
            # -unoptimized structure) -- mirrors LinkedMolOrcaFW, but with the
            # Gaussian-style mixed get_from_run_dict/get_from_str mol dict
            # that class doesn't support, so it is built directly here.
            node_working_dir = os.path.join(
                working_dir, self.dir_head, *dir_structure
            )
            link_task = ProcessMoleculeInput(
                mol=self.link_mol,
                operation_type="link_molecules",
                from_fw_spec=True,
                str_type="smi",
                db=db,
            )
            if not self.skip:
                opt_task = RunOrca(
                    gout_key=self.gout_key + "_opt",
                    db=db,
                    solvent=node_solvent,
                    tag=tag,
                    **{
                        i: j
                        for i, j in {**opt_orca_inputs, **orca_settings}.items()
                        if i in RunOrca.optional_params
                    },
                )
                link_opt_fw = Firework(
                    [link_task, opt_task],
                    name=f"{self.mol_name}_optimization",
                    spec={
                        "tag": tag,
                        "_launch_dir": os.path.join(node_working_dir, "Optimization"),
                    },
                )
                freq_fw = OrcaFW(
                    prev_calc_key=self.gout_key + "_opt",
                    gaussian_input_params=freq_orca_inputs,
                    db=db,
                    name=f"{self.mol_name}_frequency",
                    parents=link_opt_fw,
                    working_dir=os.path.join(node_working_dir, "Frequency"),
                    gout_key=self.gout_key,
                    tag=tag,
                    solvent=node_solvent,
                    **orca_settings,
                )
                self.fireworks = [link_opt_fw, freq_fw]
            else:
                freq_task = RunOrca(
                    gout_key=self.gout_key,
                    db=db,
                    solvent=node_solvent,
                    tag=tag,
                    **{
                        i: j
                        for i, j in {**freq_orca_inputs, **orca_settings}.items()
                        if i in RunOrca.optional_params
                    },
                )
                link_freq_fw = Firework(
                    [link_task, freq_task],
                    name=f"{self.mol_name}_frequency",
                    spec={
                        "tag": tag,
                        "_launch_dir": os.path.join(node_working_dir, "Frequency"),
                    },
                )
                self.fireworks = [link_freq_fw]
        else:
            mol_kwargs = (
                {"mol": self.mol} if self.parent is None else {"prev_calc_key": self.parent.gout_key}
            )
            _, label, local_fws = common_fw(
                working_dir=working_dir,
                opt_gaussian_inputs=opt_orca_inputs,
                freq_gaussian_inputs=freq_orca_inputs,
                gout_key=self.gout_key,
                db=db,
                mol_name=self.mol_name,
                dir_head=self.dir_head,
                dir_structure=dir_structure,
                skips=self.skip,
                tag=tag,
                solvent=node_solvent,
                **mol_kwargs,
                **orca_settings,
            )
            if self.parent is None and self.dir_head is None:
                self.dir_head = label
            self.fireworks = local_fws
        return

    def branch(
        self,
        branching_states,
        branching_phases,
        num_of_electrons,
        branch_cation_from_anion,
        h_index,
        vertical,
    ):
        """Generate the children nodes of the current node in the tree representing the IP/EA workflow."""
        if self.state == "cation":
            branching_states = [i for i in branching_states if i != "anion"]
        if self.state == "anion":
            if not branch_cation_from_anion:
                branching_states = [i for i in branching_states if i != "cation"]
        for state in branching_states:
            if vertical:
                skips = None if state == "reference" else ["opt"]
            else:
                skips = None
            self.children_nodes.append(
                Node(
                    state,
                    self.phase,
                    num_of_electrons,
                    parent=self,
                    branch_cation_from_anion=branch_cation_from_anion,
                    h_index=h_index,
                    skips=skips,
                )
            )
        if self.phase in branching_phases:
            branching_phases.remove(self.phase)
        for phase in branching_phases:
            if vertical:
                skips = None if self.state == "reference" else ["opt"]
            else:
                skips = None
            self.children_nodes.append(
                Node(
                    self.state,
                    phase,
                    0,
                    parent=self,
                    branch_cation_from_anion=branch_cation_from_anion,
                    skips=skips,
                )
            )
        return self.children_nodes


def _validate_ip_ea_inputs(ref_charge, opt_orca_inputs, states, phases, pcet, h_index, num_electrons):
    """Resolve states/phases defaults and validate the ip_ea input combination; raises on any inconsistency."""
    if ref_charge != opt_orca_inputs.get("charge", ref_charge):
        raise Exception(
            "The provided reference charge is not consistent with "
            "the one found in the orca input parameters."
        )

    if states is None:
        states = ["cation", "anion"]
    if phases is None:
        phases = ["gas", "solution"]

    if not states:
        raise ValueError("states list is empty")
    for state in states:
        if state.lower() not in ["cation", "anion"]:
            raise ValueError(
                "The provided states are not supported. Supported"
                " ones are reference, cation, and/or anion."
            )
    if not phases:
        raise ValueError("phases list is empty")
    for phase in phases:
        if phase.lower() not in ["gas", "solution"]:
            raise ValueError(
                "The provided phases are not supported. Supported"
                " ones are gas and/or solution."
            )

    if pcet:
        if h_index is None:
            raise ValueError(
                "index at which to attach hydrogen atom should be provided as "
                "input"
            )
        if len(h_index) != num_electrons:
            raise ValueError(
                "number of indices at which to attach hydrogen atoms should be "
                "consistent with number of transfer steps"
            )

    return states, phases


def _normalize_electrode_potentials(electrode_potentials, working_dir):
    """Lower-case the electrode_potentials keys and resolve each entry's "ref" bibtex key to its full citation."""
    if not electrode_potentials:
        return electrode_potentials

    electrode_potentials = {
        k.lower(): {i.lower(): j for i, j in v.items()}
        if isinstance(v, dict)
        else v
        for k, v in electrode_potentials.items()
    }
    for k, v in electrode_potentials.items():
        if type(v) != dict or "potential" and "ref" not in v:
            raise KeyError(
                "Standard electrode potential dict should "
                "contain potential and ref keys."
            )
        electrode_potentials[k]["ref"] = bibtex_parser(v["ref"], working_dir)
    return electrode_potentials


def _build_orca_settings(orca_cmd, num_cores, memory):
    """Collect the non-default engine-level ORCA settings into the kwargs dict RunOrca/OrcaFW expect."""
    orca_settings = {}
    if orca_cmd:
        orca_settings["orca_cmd"] = orca_cmd
    if num_cores:
        orca_settings["num_cores"] = num_cores
    if memory:
        orca_settings["memory"] = memory
    return orca_settings


def _branch_from_node(current_node, states, phases, pcet, h_index, num_electrons, single_step, vertical):
    """Compute the branching states/phases for current_node's children and create them, per the tree-growth rules."""
    if single_step:
        addition_electrons = num_electrons
    else:
        addition_electrons = 1

    if not (
        abs(current_node.added_e) <= num_electrons
        or abs(current_node.added_h) <= num_electrons
    ):
        return []

    if "gas" in phases:
        if current_node.phase == "solution":
            branching_phases = []
            branching_states = []
        else:
            branching_phases = ["solution"]
            branching_states = deepcopy(states)
    else:
        branching_phases = []
        branching_states = deepcopy(states)

    if abs(current_node.added_e) == num_electrons:
        if pcet and current_node.added_h < num_electrons:
            branching_states = [i for i in branching_states if i == "cation"]
        else:
            branching_states = []
    if abs(current_node.added_h) == num_electrons:
        branching_states = []

    branching_phases = [i for i in branching_phases if i in phases]
    return current_node.branch(
        branching_states,
        branching_phases,
        addition_electrons,
        pcet,
        h_index,
        vertical,
    )


def get_ip_ea(
    mol_operation_type,
    mol,
    ref_charge,
    single_step=False,
    vertical=False,
    pcet=False,
    h_index=None,
    num_electrons=1,
    opt_orca_inputs=None,
    freq_orca_inputs=None,
    solvent_gaussian_inputs=None,
    solvent_properties=None,
    states=None,
    phases=None,
    electrode_potentials=None,
    gibbs_elec=-0.001378786,
    gibbs_h=-0.41816,
    db=None,
    name="ip_ea_calculation",
    working_dir=None,
    orca_cmd=None,
    num_cores=None,
    memory=None,
    ref_skips=None,
    **kwargs,
):
    """Define a workflow for calculating the ionization potential (IP) and electron affinity (EA) in eV, using ORCA.

    Supports the same methods as
    ``mispr.gaussian.workflows.base.ip_ea.get_ip_ea``:

    * **Direct electron transfer**
    * **Vertical calculation of IP and EA**
    * **Adiabatic IP/EA**
    * **Sequential PCET**

    Uses the same tree structure as the Gaussian workflow to dynamically
    define the dependencies of the Fireworks; see that function's docstring
    for the physical meaning of each argument. Differences here are purely
    ORCA-backend plumbing:

    Args:
        mol_operation_type (str): The type of molecule operation; see
            ``process_mol`` in ``mispr/gaussian/utilities/mol.py``.
        mol (Molecule, str, dict): Source of the molecule to be processed;
            should match ``mol_operation_type``.
        ref_charge (int): The initial charge on the molecule.
        single_step (bool, optional): Defaults to ``False``.
        vertical (bool, optional): Defaults to ``False``.
        pcet (bool, optional): Defaults to ``False``.
        h_index (list, optional): Site indices for PCET hydrogen attachment.
        num_electrons (int, optional): Defaults to 1.
        opt_orca_inputs (dict, optional): ORCA input parameters for the
            optimization step, in the same shape as ``OrcaFW``'s
            ``gaussian_input_params`` (e.g. {"functional": "b3lyp",
            "basis_set": "6-31G(d)", "route_parameters": {"Opt": None}});
            defaults to B3LYP/6-31G(d).
        freq_orca_inputs (dict, optional): Same, for the frequency step;
            defaults to B3LYP/6-31G(d).
        solvent_gaussian_inputs (str, optional): Gaussian-style implicit
            solvent string (e.g. "(Solvent=Water)"), kept in this format for
            consistency with the Gaussian workflow's db metadata and
            translated internally into ORCA's CPCM keyword; defaults to
            "(PCM, Solvent=Water)" if "solution" is in ``phases``.
        solvent_properties (dict, optional): Recorded in the final db document
            only; ORCA's CPCM model has no equivalent input knob for this.
        states (list, optional): Defaults to ["cation", "anion"].
        phases (list, optional): Defaults to ["gas", "solution"].
        electrode_potentials (dict, optional): See the Gaussian workflow.
        gibbs_elec (float, optional): Electron Gibbs free energy in Hartree.
        gibbs_h (float, optional): Hydrogen Gibbs free energy in Hartree.
        db (str or dict, optional): Database credentials.
        name (str, optional): Name of the workflow.
        working_dir (str, optional): Working directory; defaults to cwd.
        orca_cmd (str, optional): Path to the ORCA executable; falls back to
            the ORCA_CMD environment variable, then to "orca" on PATH.
        num_cores (int, optional): Number of parallel ORCA processes per
            calculation; defaults to 1.
        memory (int, optional): Memory per core in MB; defaults to 4000.
        ref_skips (list, optional): Jobs to skip for the reference (neutral)
            state only; e.g. ["opt"].
        kwargs (keyword arguments): Additional kwargs forwarded to
            ``IPEAtoDB``/``Workflow`` (e.g. ``tag``).

    Returns:
        Workflow

    """
    fws = []
    fireworks_dict = {}
    links_dict = {}
    parents_dict = {}
    working_dir = working_dir or os.getcwd()
    mol = recursive_relative_to_absolute_path(mol, working_dir)

    opt_orca_inputs = opt_orca_inputs or {
        "functional": "b3lyp",
        "basis_set": "6-31G(d)",
        "route_parameters": {"Opt": None},
    }
    freq_orca_inputs = freq_orca_inputs or {
        "functional": "b3lyp",
        "basis_set": "6-31G(d)",
        "route_parameters": {"Freq": None},
    }

    states, phases = _validate_ip_ea_inputs(
        ref_charge, opt_orca_inputs, states, phases, pcet, h_index, num_electrons
    )

    if "solution" in phases and not solvent_gaussian_inputs:
        solvent_gaussian_inputs = "(PCM, Solvent=Water)"
    solvent = _to_orca_solvent(solvent_gaussian_inputs)

    electrode_potentials = _normalize_electrode_potentials(
        electrode_potentials, working_dir
    )

    orca_settings = _build_orca_settings(orca_cmd, num_cores, memory)
    tag = kwargs.get("tag", "unknown")

    root_mol = process_mol(mol_operation_type, mol, db=db)
    root_node = Node(
        "reference",
        "gas" if "gas" in phases else "solution",
        0,
        mol=root_mol,
        skips=ref_skips,
        ref_charge=ref_charge,
        branch_cation_from_anion=pcet,
    )
    solved_nodes = []
    active_nodes = Queue()
    active_nodes.put(root_node)
    while not active_nodes.empty():
        current_node = active_nodes.get()
        current_node.create_fireworks(
            opt_orca_inputs,
            freq_orca_inputs,
            solvent,
            working_dir,
            db,
            pcet,
            tag,
            orca_settings,
            **kwargs,
        )
        fws += current_node.fireworks
        fireworks_dict[current_node.gout_key] = current_node.fireworks
        if current_node.parent is not None:
            links_dict[current_node.parent.fireworks[-1]] = links_dict.get(
                current_node.parent.fireworks[-1], []
            ) + [current_node.fireworks[0]]
            parents_dict[current_node.gout_key] = current_node.parent.gout_key

        children = _branch_from_node(
            current_node, states, phases, pcet, h_index, num_electrons, single_step, vertical
        )
        for child in children:
            active_nodes.put(child)
        solved_nodes.append(current_node)

    gout_keys = [i.gout_key for i in solved_nodes]
    fw_analysis = Firework(
        IPEAtoDB(
            num_electrons=num_electrons,
            states=states,
            phases=phases,
            steps="single" if single_step else "multi",
            root_node_key=root_node.gout_key,
            keys=gout_keys,
            pcet=pcet,
            vertical=vertical,
            solvent_gaussian_inputs=solvent_gaussian_inputs,
            solvent_properties=solvent_properties,
            electrode_potentials=electrode_potentials,
            gibbs_elec=gibbs_elec,
            gibbs_h=gibbs_h,
            db=db,
            **{
                i: j
                for i, j in kwargs.items()
                if i in IPEAtoDB.required_params + IPEAtoDB.optional_params
            },
        ),
        parents=fws[:],
        name="{}-{}".format(root_node.dir_head, "ip_ea_analysis"),
        spec={
            "tag": tag,
            "_launch_dir": os.path.join(working_dir, root_node.dir_head, "analysis"),
        },
    )
    fws.append(fw_analysis)

    return Workflow(
        fws,
        name="{}_{}".format(root_node.dir_head, name),
        links_dict=links_dict,
        **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS},
    )
