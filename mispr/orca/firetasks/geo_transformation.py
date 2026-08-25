"""Define ORCA-flavored geometry transformation firetasks.

The bond-breaking/fragmentation logic in
``mispr.gaussian.firetasks.geo_transformation.BreakMolecule`` (charge assignment,
finding unique fragments, ring opening, etc.) has nothing to do with which quantum
chemistry engine is used downstream -- only the dynamically generated per-fragment
optimization/frequency Fireworks are engine-specific. This module therefore
subclasses the Gaussian ``BreakMolecule`` and overrides only the one method
(``_workflow``) that decides which Fireworks to build for each fragment, pointing
it at the ORCA ``common_fw`` instead of the Gaussian one.
"""

import logging

from fireworks import Workflow
from fireworks.utilities.fw_utilities import explicit_serialize

from mispr.gaussian.utilities.metadata import get_mol_formula
from mispr.gaussian.firetasks.geo_transformation import (
    BreakMolecule as GaussianBreakMolecule,
)

__author__ = "Ruiqi Luo"
__status__ = "Development"
__date__ = "2026_7_19"
__version__ = "0.0.5"

logger = logging.getLogger(__name__)

# engine-level RunOrca settings that must survive the kwargs filtering below
# and reach the per-fragment Fireworks
ORCA_ENGINE_KWARGS = ("orca_cmd", "num_cores", "memory")


@explicit_serialize
class BreakMolecule(GaussianBreakMolecule):
    """ORCA counterpart of ``mispr.gaussian.firetasks.geo_transformation.BreakMolecule``.

    Identical fragmentation logic, but generates ORCA optimization/frequency
    Fireworks for each fragment instead of Gaussian ones.
    """

    @staticmethod
    def _workflow(
        mol,
        gout_key,
        working_dir,
        db,
        opt_gaussian_inputs,
        freq_gaussian_inputs,
        cart_coords,
        oxidation_states,
        save_to_db,
        save_to_file,
        fmt,
        update_duplicates,
        **kwargs,
    ):
        """Build the ORCA opt/freq Fireworks for one already-charged/spin-set fragment molecule.

        All arguments match the base class's abstract ``_workflow`` signature
        (see ``mispr.gaussian.firetasks.geo_transformation.BreakMolecule`` for
        what each one means) -- only the body differs here, building ORCA
        Fireworks via ``common_fw`` instead of Gaussian ones. Single-atom
        fragments have their "opt" job skipped (a lone atom has no geometry to
        optimize).

        Returns:
            Workflow

        """
        from mispr.orca.workflows.base.core import common_fw, WORKFLOW_KWARGS

        dir_structure = ["charge_{}".format(str(mol.charge))]
        mol_formula = get_mol_formula(mol)

        if len(mol) == 1:
            skips = ["opt"]
        else:
            skips = None

        job_name = "opt_freq"
        _, _, frag_fws = common_fw(
            mol=mol,
            working_dir=working_dir,
            dir_structure=dir_structure,
            db=db,
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            mol_name=mol_formula,
            skips=skips,
            gout_key=gout_key,
            **{
                i: j
                for i, j in kwargs.items()
                if i in WORKFLOW_KWARGS
                or i == "tag"
                or i in ORCA_ENGINE_KWARGS
            },
        )
        return Workflow(
            frag_fws,
            name="{}_{}".format(mol_formula, job_name),
            **{i: j for i, j in kwargs.items() if i in WORKFLOW_KWARGS},
        )
