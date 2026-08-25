===========
Keywords
===========

.. list-table::
   :widths: 25 75
   :header-rows: 1

   * - Keyword
     - Definition
   * - ``|CODES_DIR|``
     - Main directory where the virtual python environment is created and the configuration files are stored
   * - ``LaunchPad``
     - FireWorks database that controls the workflows. It stores all the tasks to be run and their status (e.g., RUNNING, WAITING, COMPLETED, etc.)
   * - ``Workflow``
     - A Directed Acyclic Graph of one or more Fireworks with dependencies between them, describing the full procedure for computing a property
   * - ``FireWork``
     - A list of FireTasks that are to be run in sequence; one node in a Workflow (e.g. "Optimization", "Frequency", "ESP")
   * - ``FireTask``
     - Computing task to be performed; the smallest unit of work in a Firework (e.g. write an input file, submit a job, parse an output file)
   * - ``FireWorker``
     - Defines where and how a job runs -- can be as simple as the local workstation hosting the LaunchPad, or a supercomputing cluster with a queueing system; configured in ``my_fworker.yaml``
   * - ``Rocket`` / ``rlaunch``
     - The process that pulls a FireWork off the LaunchPad and executes it; ``rlaunch`` runs one Rocket, ``qlaunch`` submits Rockets through a queueing system
   * - ``runs`` collection
     - The MongoDB collection where MISPR stores one document per individual calculation step; distinct from the property-specific collections (e.g. ``esp``, ``bde``) that hold the final, combined-analysis document for a workflow
   * - ``orca_cmd`` / ``ORCA_CMD``
     - How the ORCA backend locates the ``orca`` executable: the ``orca_cmd`` argument accepted by every ORCA workflow function, or the ``ORCA_CMD`` environment variable if ``orca_cmd`` is not given. Unlike Gaussian/LAMMPS, this is not read from ``config.ini``
