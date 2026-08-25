"""Define firework used to break a molecule and run its fragments with ORCA."""

import os
import logging

from fireworks import Firework

from mispr.gaussian.firetasks.geo_transformation import ProcessMoleculeInput
from mispr.orca.firetasks.geo_transformation import BreakMolecule

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames


class BreakMolFW(Firework):
    """ORCA counterpart of ``mispr.gaussian.fireworks.break_mol.BreakMolFW``.

    Process a molecule input, break it into unique fragments, and generate a
    set of ORCA optimization and frequency calculations for each fragment
    (optional).

    ``ProcessMoleculeInput`` (molecule-format handling) is reused unmodified
    from the Gaussian firetasks since it has no dependency on the QM engine
    used downstream; only ``BreakMolecule`` (which decides what Fireworks to
    spawn per fragment) is ORCA-specific.
    """

    def __init__(
        self,
        mol,
        mol_operation_type="get_from_mol",
        bonds=None,
        open_rings=False,
        ref_charge=0,
        fragment_charges=None,
        calc_frags=True,
        db=None,
        name="break_mol",
        parents=None,
        working_dir=None,
        tag="unknown",
        **kwargs,
    ):
        """See ``mispr.gaussian.fireworks.break_mol.BreakMolFW`` for the full list of arguments; identical signature and semantics.

        Engine-level ORCA settings (``orca_cmd``, ``num_cores``, ``memory``)
        passed via kwargs are forwarded to the per-fragment Fireworks.
        """
        t = []
        working_dir = working_dir or os.getcwd()
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        t.append(
            ProcessMoleculeInput(
                mol=mol,
                operation_type=mol_operation_type,
                db=db,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in ProcessMoleculeInput.required_params
                    + ProcessMoleculeInput.optional_params
                },
            )
        )

        t.append(
            BreakMolecule(
                bonds=bonds,
                open_rings=open_rings,
                ref_charge=ref_charge,
                fragment_charges=fragment_charges,
                calc_frags=calc_frags,
                db=db,
                additional_kwargs=kwargs,
                **{
                    i: j
                    for i, j in kwargs.items()
                    if i
                    in BreakMolecule.required_params + BreakMolecule.optional_params
                },
            )
        )

        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(BreakMolFW, self).__init__(
            t,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS},
        )
