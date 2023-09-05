# coding: utf-8


# Defines firework used to break a molecule and run its fragments.

import os
import logging

from fireworks import Firework

from mispr.gaussian.firetasks.geo_transformation import (
    BreakMolecule,
    ProcessMoleculeInput,
)

__author__ = "Rasha Atwi"
__maintainer__ = "Rasha Atwi"
__email__ = "rasha.atwi@stonybrook.edu"
__status__ = "Development"
__date__ = "Jan 2021"
__version__ = "0.0.4"

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames


class BreakMolFW(Firework):
    """
    Processes a molecule input, breaks it into unique fragments, and
    generates a set of optimization and frequency calculations for
    each fragment (optional).
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
        **kwargs
    ):
        """
        Args:
            mol (Molecule, GaussianOutput, str, dict): source of the
                molecule to be processed. Should match the mol_operation_type
            mol_operation_type (str): the type of molecule operation.
                See process_mol defined in mispr/gaussian/utilities/mol.py
                for supported operations; defaults to "get_from_mol"
            bonds (list): list of tuples of the bonds to break; e.g.
                [(0, 1), (1, 2)] will break the bonds between atoms 0
                 and 1 and between atoms 1 and 2; if none is specified,
                 will attempt to break all bonds
            open_rings (bool): whether to open rings; if set to True, will
                perform local optimization to get a good initial guess for
                the structure; defaults to False
            ref_charge (int): charge on the principle molecule;
                defaults to 0
            fragment_charges (list): list of charges to assign to
                the fragments in addition to the ones already assigned;
                refer to mispr.gaussian.firetasks.geo_transformation.BreakMolecule
                for more details
            calc_frags (bool): whether to create optimization and frequency
                Fireworks for the generated fragments; defaults to True
            db (str or dict): database credentials
            name (str): name of the Firework; defaults to "break_mol"
            parents (Firework or [Firework]): list of parent FWs this FW
                depends on
            working_dir (str): working directory for the calculation;
                will use the current working directory if not specified
            tag (str): tag for the calculation; the provided tag will be
                stored in the db documents for easy retrieval; defaults
                to "unknown"
            **kwargs: other kwargs that are passed to:
                1. Firework.__init__.
                2. mispr.gaussian.firetasks.geo_transformation.ProcessMoleculeInput
                3. mispr.gaussian.firetasks.geo_transformation.BreakMolecule
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
                }
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
                }
            )
        )

        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(BreakMolFW, self).__init__(
            t,
            parents=parents,
            name=name,
            spec=spec,
            **{i: j for i, j in kwargs.items() if i in FIREWORK_KWARGS}
        )
