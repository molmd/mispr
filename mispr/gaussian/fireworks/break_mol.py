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
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

FIREWORK_KWARGS = Firework.__init__.__code__.co_varnames


class BreakMolFW(Firework):
    def __init__(self,
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
                 **kwargs):
        t = []
        working_dir = working_dir or os.getcwd()
        if not os.path.exists(working_dir):
            os.makedirs(working_dir)

        t.append(ProcessMoleculeInput(mol=mol,
                                      operation_type=mol_operation_type,
                                      db=db,
                                      **{i: j for i, j in kwargs.items() if i in
                                         ProcessMoleculeInput.required_params +
                                         ProcessMoleculeInput.optional_params}
                                      )
                 )

        t.append(BreakMolecule(bonds=bonds,
                               open_rings=open_rings,
                               ref_charge=ref_charge,
                               fragment_charges=fragment_charges,
                               calc_frags=calc_frags,
                               db=db,
                               additional_kwargs=kwargs,
                               **{i: j for i, j in kwargs.items() if i in
                                  BreakMolecule.required_params +
                                  BreakMolecule.optional_params}
                               )
                 )

        spec = kwargs.pop("spec", {})
        spec.update({"tag": tag, "_launch_dir": working_dir})
        super(BreakMolFW, self).__init__(t,
                                         parents=parents,
                                         name=name,
                                         spec=spec,
                                         **{i: j for i, j in kwargs.items()
                                            if i in FIREWORK_KWARGS})
