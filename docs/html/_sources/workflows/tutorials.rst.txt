===================
Workflow Tutorials
===================

This page is under construction.

Running an ESP workflow
------------------------------
This workflow calculates the electrostatic partial charges that are fit to the electrostatic potential at points
selected according to the Merz-Singh-Kollman scheme, but other schemes supported by Gaussian can be used.
The ESP workflow performs the following steps:


.. mermaid::

    flowchart TD
        A{Input Structure} -->|Preprocessing| B(Geometry Optimization)
        B --> C(Frequency Calculations)
        C --> D[ESP Calculation]
        D --> E(Analysis)





Running an MD workflow
------------------------------


Running a hybrid workflow
------------------------------
