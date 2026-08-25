====================
Supported Workflows
====================
Some of the workflows available as of July 2022 are:

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

One can customize any of the above workflows or create their own by reusing
the building blocks provided by MISPR. The above preset workflows are in
``mispr/gaussian/workflows/base``, ``mispr/lammps/workflows/base``, and
``mispr/hybrid/workflows``.

.. note::
    Other types of force field parameters can be provided as
    inputs to the MD workflow, thereby skipping the force field
    generation step.