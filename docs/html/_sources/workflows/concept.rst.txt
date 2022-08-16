==========
Workflow Concept
==========
A scientific workflow in MISPR provides a complete description of the
procedure leading to the final data used to predict the desired property
of a given molecule or system. It consists of multiple steps ranging
from the initial setup of a molecule or system of molecules to a
sequence of calculations with dependencies and optional automated
post-processing of parsed data to derive properties of interest.
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
An example of the structure of a workflow in MISPR is shown below:

.. figure:: ../static/workflow.png

Once a Workflow object is created, the user can use the FireWorks package
to execute the calculations on various computing resources. The goal of
MISPR infrastructure is to provide preset workflows for
computing properties relevant to the molecular science community and to
simplify the process of creating new workflows by using the implemented
FireWorks and Firetasks in MISPR as building blocks for custom workflows.
