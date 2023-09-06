===================
Workflow Tutorials
===================

This page is under construction.

Running an ESP workflow
------------------------------
This workflow calculates the partial charges on atoms of a molecule. The charges are fit to the electrostatic potential at
points selected according to the Merz-Singh-Kollman scheme, but other schemes supported by Gaussian can be used as well.

**The ESP workflow performs the following steps:**


.. mermaid::

    graph TD
        A[(Input Structure)] -->|Preprocessing| DFT
        DFT -->| | B[Geometry Optimization]
        B -->| | C[Frequency Calculation]
        C -->| | D[ESP Calculation]
        D -->|Postprocessing| E[(Output)]

        subgraph DFT
        B[Geometry Optimization]
        C[Frequency Calculation]
        D[ESP Calculation]
        end

        style A fill:#EBEBEB, stroke:#BB2528, arrowColor:#A9A9A9
        style DFT fill:#DDEEFF,stroke:#DDEEFF,font-weight:bold
        style B fill:#fff,stroke-dasharray: 5, 5, stroke:#BB2528
        style C fill:#fff,stroke-dasharray: 5, 5, stroke:#BB2528
        style D fill:#fff,stroke:#BB2528
        style E fill:#EBEBEB,stroke:#BB2528


In the following example, we will run the ESP workflow on a water molecule.
The input structure is provided in the file `water.xyz`:

.. code-block:: python
    :linenos:

    import os

    from fireworks import LaunchPad

    from mispr.gaussian.workflows.base.esp import get_esp_charges

    lpad = LaunchPad.auto_load()
    wf, _ = get_esp_charges(
        mol_operation_type="get_from_pubchem",
        mol="water",
        format_chk=True,
        save_to_db=True,
        save_to_file=True,
        additional_prop_doc_fields={"solvent": "water"},
        tag="mispr_tutorial",
    )
    lpad.add_wf(wf)

.. admonition:: Note
   :class: toggle

   Here, we explain the specific lines of code:

   - Line 3: Imports the necessary modules.
   - Line 5: Loads the LaunchPad.
   - Lines 7-18: Generates an ESP charges workflow for water.
   - Line 19: Adds the workflow to the LaunchPad.






Running an MD workflow
------------------------------


Running a hybrid workflow
------------------------------
