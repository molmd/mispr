===================
Workflow Tutorials
===================

This page is under construction.

Running an ESP workflow
------------------------------
This workflow calculates the partial charges on atoms of a molecule. The charges are fit to the electrostatic potential at
points selected according to the Merz-Singh-Kollman scheme, but other schemes supported by Gaussian can be used as well.

This workflow calculates the partial charges on atoms of a molecule. The charges are fit to the electrostatic potential at
points selected according to the Merz-Singh-Kollman scheme, but other schemes supported by Gaussian can be used as well.

**The ESP workflow performs the following steps:**


.. mermaid::

    %%{
    init: {
        'theme': 'base',
        'themeVariables': {
        'primaryTextColor': 'black',
        'lineColor': 'lightgrey',
        'secondaryColor': 'pink',
        'tertiaryColor': 'lightgrey'
        }
    }
    }%%

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

        style A fill:#EBEBEB,stroke:#BB2528
        style DFT fill:#DDEEFF,stroke:#DDEEFF,font-weight:bold
        style B fill:#fff,stroke-dasharray: 5, 5, stroke:#BB2528
        style C fill:#fff,stroke-dasharray: 5, 5, stroke:#BB2528
        style D fill:#fff,stroke:#BB2528
        style E fill:#EBEBEB,stroke:#BB2528


In the following example, we will run the ESP workflow on a monoglyme molecule.

.. code-block:: python
    :linenos:

    import os

    from fireworks import LaunchPad

    from mispr.gaussian.workflows.base.esp import get_esp_charges

    lpad = LaunchPad.auto_load()
    wf, _ = get_esp_charges(
        mol_operation_type="get_from_pubchem", # (1)!
        mol="monoglyme",
        format_chk=True,
        save_to_db=True,
        save_to_file=True,
        additional_prop_doc_fields={"solvent": "water"},
        tag="mispr_tutorial",
    )
    lpad.add_wf(wf)

.. code-annotations::
    1.
        :code:`mol_operation_type` refers to the operation to be performed on the input to process the molecule.

        In this example, we are requesting to directly retrieve the molecule from PubChem by providing a
        common name for the molecule to be used as query criteria for searching the PubChem database via
        the :code:`mol` input argument. For a list of supported :code:`mol_operation_type` and the corresponding
        :code:`mol`, refer to :doc:`installation guide <../mispr.gaussian.utilities.mol.process_mol>`.


Download :download:`esp_tutorial.py <../_downloads/esp_tutorial.py>`.


Running an MD workflow
------------------------------


Running a hybrid workflow
------------------------------
