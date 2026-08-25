====================
Supported Workflows
====================
The available preset workflows are:

* DFT:
   * Electrostatic partial charges (ESP)
   * NMR shifts
   * Redox potentials
   * Binding energies
   * Bond dissociation energies
* MD:
   * Initial configuration building, generation of `GAFF <http://ambermd.org>`_
     or `OPLS <http://zarbi.chem.yale.edu/oplsaam.html>`_ parameters,
     running of MD simulations
   * Analysis of output and trajectory files (e.g. RDF, coordination
     number, diffusion coefficients, etc.)
* Hybrid:
   * Core workflow for optimizing the individual structure of the
     mixture of components, generating their ESP charges, and using
     them in MD simulations
   * NMR: deriving NMR chemicals for stable solvation structures
     extracted from MD simulations

The DFT workflows run with `Gaussian <https://gaussian.com>`_ by default. As
of version 0.0.5, the ESP, bond dissociation energy, and binding energy
workflows are additionally available with an
`ORCA <https://www.faccts.de/orca/>`_ backend
(``mispr/orca/workflows/base``); the workflow functions keep the same names
and arguments across both engines, so switching engine is an import-path
change, not a code change. See
:doc:`Workflow Tutorials <tutorials>` for per-backend setup and examples.

One can customize any of the above workflows or create their own by reusing
the building blocks provided by MISPR. The preset workflows are in
``mispr/gaussian/workflows/base``, ``mispr/orca/workflows/base``,
``mispr/lammps/workflows/base``, and ``mispr/hybrid/workflows``.

.. note::
    Other types of force field parameters can be provided as
    inputs to the MD workflow, thereby skipping the force field
    generation step.