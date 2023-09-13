===================
Workflow Tutorials
===================

This page is under construction.

Running an ESP workflow
------------------------------
The ESP workflow calculates the partial charges on atoms of a molecule. The charges are
fit to the electrostatic potential at points selected according to the Merz-Singh-Kollman
scheme, but other schemes supported by Gaussian can be used as well.

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

.. note::
    The geometry optimization and frequency calculation steps (marked with a dashed
    border in the above diagram) are optional. If the input structure is already
    optimized, the workflow will skip these steps.


In the following example, we will run the ESP workflow on a monoglyme molecule.

.. code-block:: python
    :linenos:

    from fireworks import LaunchPad

    from mispr.gaussian.workflows.base.esp import get_esp_charges

    lpad = LaunchPad.auto_load()

    wf, _ = get_esp_charges(
        mol_operation_type="get_from_pubchem", # (1)!
        mol="monoglyme",
        format_chk=True,
        save_to_db=True,
        save_to_file=True,
        additional_prop_doc_fields={"name": "monoglyme"},
        tag="mispr_tutorial",
    )
    lpad.add_wf(wf) # (2)!

.. code-annotations::
    1.
        :code:`mol_operation_type` refers to the operation to be performed on the input to process the molecule.

        In this example, we are requesting to directly retrieve the molecule from PubChem by providing a
        common name for the molecule to be used as query criteria for searching the PubChem database via
        the :code:`mol` input argument. For a list of supported :code:`mol_operation_type` and the corresponding
        :code:`mol`, refer to :meth:`mispr.gaussian.utilities.mol.process_mol`.

    2. Adds the workflow to the launchpad.


Download :download:`esp_tutorial.py <../_downloads/esp_tutorial.py>`.

Run the script using the following command:

.. code-block:: bash

    python esp_tutorial.py

And then launch the job through the queueing system using the following command:

.. code-block:: bash

    qlaunch rapidfire # (1)!

.. code-annotations::
    1.
        This command can submit a large number of jobs at once
        or maintain a certain number of jobs in the queue.

The workflow will run and create a directory named :code:`C4H10O2` in the current working
directory. The directory will contain the following subdirectories:

.. code-block:: bash

    C4H10O2
    ├── Optimization
    ├── Frequency
    ├── ESP
    ├── analysis

Inside the :code:`Optimization`, :code:`Frequency`, and :code:`ESP` subdirectories, you
will find the Gaussian input and output files for the corresponding step. Inside the
:code:`Optimization` subdirectory, you will also find a "convergence.png" figure that
shows the forces and displacement convergence during the course of the optimization.

.. figure:: ../_static/convergence.png

The :code:`analysis` subdirectory contains the results of the workflow in the form of a
:code:`esp.json` file. You can read the content of the :code:`esp.json` file using the
following commands:

.. code-block:: python
    :linenos:

    import json

    with open("C4H10O2/analysis/esp.json", "r") as f:
        esp = json.load(f)

    print(esp["esp"])

This will output the partial charges on the atoms of the molecule:

.. code-block:: python

    {
    "1": ["O", -0.374646],
    "2": ["O", -0.373831],
    "3": ["C", 0.132166],
    "4": ["C", 0.132716],
    "5": ["C", 0.034284],
    "6": ["C", 0.031733],
    "7": ["H", 0.033853],
    "8": ["H", 0.034024],
    "9": ["H", 0.034218],
    "10": ["H", 0.034388],
    "11": ["H", 0.070724],
    "12": ["H", 0.03474],
    "13": ["H", 0.03438],
    "14": ["H", 0.034621],
    "15": ["H", 0.071656],
    "16": ["H", 0.034974],
    }

Running a BDE workflow
------------------------------


Running an MD workflow
------------------------------


Running a hybrid workflow
------------------------------
