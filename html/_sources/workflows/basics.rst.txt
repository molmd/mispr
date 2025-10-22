====================
Workflow Basics
====================
A scientific workflow in MISPR provides a complete description of the
procedure leading to the final data used to predict the desired property
of a given molecule or system. It consists of multiple steps ranging
from the initial setup of a molecule or system of molecules to a
sequence of calculations with dependencies and optional automated
post-processing of parsed data to derive properties of interest.

.. note::
   The workflow model we use to encode DFT and MD recipes in MISPR is
   defined by the FireWorks workflow software.

A workflow in FireWorks is modeled as a Directed Acyclic Graph
representing the chain of relationships between
computational operations. A workflow consists of one or more Fireworks
(jobs) with dependencies. The workflow contains information
about the links between Fireworks to execute them in the correct order.
Each Firework consists of one or more Firetasks that run sequentially.
A Firetask is an atomic computing job that can call shell scripts,
transfer files, write/delete files, or execute other Python functions.
An example of the structure of a DFT workflow in MISPR is shown below:

.. figure:: ../_static/workflow.png

Once a Workflow object is created, the user can use the FireWorks package
to execute the calculations on various computing resources. The goal of
MISPR infrastructure is to provide preset workflows for
computing properties relevant to the molecular science community and to
simplify the process of creating new workflows by using the implemented
FireWorks and Firetasks in MISPR as building blocks for custom workflows.

At the end of each workflow in MISPR, an analysis FireTask is performed
to analyze the results and generate a report. The report is in the form
of a JSON file and/or MongoDB document. It contains all the input parameters
used in the calculations, the output data, general information about the
calculation like the software version used (Gaussian, LAMMPS, MISPR, etc.),
the wall time the full run took, and chemical metadata about the molecule
or system of molecules (e.g. SMILES, InChI, molecular formula, etc.).

In general, each property predicted by MISPR workflows is the result of
multiple Gaussian or LAMMPS calculations, and the predicted property is
represented by a single file/document summarizing data and "raw" information
collected from different calculation steps. The MongoDB document
corresponding to a predicted property is stored in a MongoDB collection
named after the property. For example, bond dissociation energies are
stored in a ``bde`` collection in the database while electrostatic
partial charges are saved in an ``esp`` collection and so on. Some of
the analysis FireTasks also include optional plotting of the results.
Besides the final summary file/document, MISPR stores data from the
intermediate calculation steps into a collection called ``runs`` in the
database.

The following diagram summarizes the process in MISPR workflows to generate
the analysis files/documents:

.. figure:: ../_static/analysis.png

.. note::
   The above diagram shows one example of the structure of a workflow
   where the Fireworks are executed sequentially. Some workflows contain
   parallel Fireworks. However, the analysis Firework
   is always the last Firework in all the workflows in MISPR.