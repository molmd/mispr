.. MISPR documentation master file, created by
   sphinx-quickstart on Sat Jul 30 18:04:52 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
.. title:: mispr (Materials Science Workflows)
.. image:: _static/logo.png
   :alt: mispr logo

=======
Overview
=======
MISPR (Materials Informatics for Structure-Property Relationships) is a
high-throughput computational infrastructure aimed at guiding and
accelerating materials discovery, optimization, and deployment for
liquid solutions by seamlessly integrating density functional theory
(DFT) with classical molecular dynamics (MD) techniques. MISPR is
motivated by the Materials Genome Initiative (MGI) principles and is
built on top of open-source Python packages developed for the Materials
Project such as `pymatgen <https://pymatgen.org>`_,
`FireWorks <https://materialsproject.github.io/fireworks/>`_ ,
and `custodian <https://materialsproject.github.io/custodian/>`_.

Features of MISPR include:

* Automation of DFT and MD simulations and all their underlying tasks
  from file management and job submission to supercomputing resources,
  to output parsing and data analytics; a task that can be done to a
  single molecule/system or to a large number of systems in parallel

* Creation of computational databases of force field parameters and DFT
  and MD derived properties of molecular systems for establishing
  structure-property relations and maintaining data provenance and
  reproducibility

* Detection of the inevitable errors that occur during the simulations
  and their on-the-fly correction based on template responses that have
  been designed relying on human intuition coupled with extensive
  experience to significantly improve the success rate of high-throughput
  simulations while eliminating human intervention

* Support for flexible and well-tested DFT workflows that compute various
  properties of individual molecular species or complexes such as bond
  dissociation energy, binding energy, redox potential, and nuclear
  magnetic resonance (NMR) tensors

* Derivation of many molecular ensemble properties such as radial
  distribution functions, diffusion coefficients, viscosity, and
  conductivity of liquid solutions, which are critical to understanding
  complex inter- and intra-atomic interactions controlling the performance
  of solutions within various chemistry, biology, and materials science
  applications

* Seamless integration of DFT and MD simulations through hybrid
  workflows that enable force field generation and information flow
  between the two length scales to allow exploring wide chemical and
  parameter spaces (e.g., temperature, pressure, concentration, etc.),
  a task that can be infeasible experimentally and challenging using
  manual calculations

* Automatic extraction of hundreds of thousands of solvation structures
  from MD ensembles and their use in DFT workflows to accurately represent
  the electronic environment, which is crucial to derive reliable energetics
  and other properties such as NMR chemical shifts and redox potentials
  and match them to experimental data

**Note**: MISPR is primarily built to work with `Gaussian <https://gaussian.com>`_
electronic structure software for DFT calculation and `LAMMPS <https://www.lammps.org/#gsc.tab=0>`_
open-source software for MD simulations.

=======
Workflows
=======
Some of the workflows available as of July 2022 are:

* DFT:
   * Electrostatic partial charges (ESP)
   * NMR shifts
   * Redox potentials
   * Binding energies
   * Bond dissociation energies
* MD:
   * Initial configuration building, generation of `GAFF <http://ambermd.org>`_ parameters,
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

**Note**: An interface for the automatic generation of
`OPLS <http://zarbi.chem.yale.edu/oplsaam.html>`_ parameters
is under development.

**Note**: Other types of force field parameters can be provided as
inputs to the MD workflow, thereby skipping the force field generation step.

=======
Getting started
=======

=======
Installation
=======

=======
What's new?
=======
Track changes to atomate through the :doc:`changelog` and the GitHub
commit log for a record of changes.

=======
Citing MISPR
=======
If you find MISPR useful in your research, please consider citing the following paper:

**Paper 1 (MISPR):**
.. code-block:: bibtex
   @article{atwi2022mispr,
     title={MISPR: An automated infrastructure for high-throughput DFT and MD simulations},
     author={Atwi, Rasha and Bliss, Matthew and Makeev, Maxim and Rajput, Nav Nidhi},
     year={2022}
   }
Download as :download:`BibTeX <_static/mispr_citation.bib>`

**Paper 2 (Hybrid NMR Workflow):**
.. code-block:: bibtex
   @article{atwi2022automated,
     title={An automated framework for high-throughput predictions of NMR chemical shifts within liquid solutions},
     author={Atwi, Rasha and Chen, Ying and Han, Kee Sung and Mueller, Karl T and Murugesan, Vijayakumar and Rajput, Nav Nidhi},
     journal={Nature Computational Science},
     volume={2},
     number={2},
     pages={112--122},
     year={2022},
     publisher={Nature Publishing Group}
   }
Download as :download:`BibTeX <_static/nmr_citation.bib>`

=======
Contributing / Contact / Support
=======

=======
License
=======


.. toctree::
   :maxdepth: 2
   :caption: Documentation:

   introduction
   API Docs<modules>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
