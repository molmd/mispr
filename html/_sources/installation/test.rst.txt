===============================
Running a Test Workflow
===============================
After installing MISPR and its dependencies and setting up the configuration
files, it is important to make sure everything is working correctly.
Therefore, we will now run a very simple test workflow that optimizes the
structure of a molecule. Later in :doc:`Workflow Tutorials <../workflows/tutorials>`,
we will show how to run more complex workflows.

.. note::
    ``lpad`` and ``qlaunch`` that are used in this tutorial are part of
    FireWorks software. You can refer to FireWorks documentation if you
    need additional information.

Initialize the database
------------------------------
Initialize the database by running the following command::

    lpad reset



.. warning::
    This command should only be executed one time when you are first
    initializing the database set-up. If you reset your LaunchPad at a
    later time, you will erase all existing entries in your FireWorks
    database, which includes your fireworks, workflows, and launches
    collections.

.. note::
    Your Python environment where FireWorks is installed must be active
    before you run this command.

Running the above command will return something like this::

    Are you sure? This will RESET 0 workflows and all data. (Y/N)y
    2022-08-15 17:04:42,224 INFO Performing db tune-up
    2022-08-15 17:04:42,683 INFO LaunchPad was RESET.

Add a workflow
------------------------------
The next step is to add a workflow to the database. We will run a
workflow that optimizes the geometry of a monoglyme molecule starting
from its xyz file. Note that you need to have the **monoglyme.xyz** file in
your working directory. You will need to run the following Python code
by creating a file called ``optimize_geometry.py``:

.. code-block:: python

    from mispr.gaussian.fireworks.core import CalcFromMolFW
    from fireworks import LaunchPad, Workflow

    lpad = LaunchPad.auto_load()
    wf = Workflow([CalcFromMolFW("monoglyme.xyz", "get_from_file",
                                  gaussian_input_params={"route_parameters": {"opt": None}},
                                  save_to_file=True, save_to_db=True)])
    lpad.add_wf(wf)

and then running the following command in terminal::

    python optimize_geometry.py


This will add a structure optimization workflow to the database.

Verify the workflow
------------------------------
To check the status of this workflow in the database, run the following
command in terminal::

    lpad get_fws -s READY

It will return something like this::

    {
        "fw_id": 1,
        "created_on": "2022-08-16T20:32:54.554404",
        "updated_on": "2022-08-16T20:32:54.554716",
        "state": "READY",
        "name": "calc_from_mol"
    }

Alternatively, you can query your ``fireworks`` collection in the MongoDB
database directly or start FireWorks' `LaunchPad <https://materialsproject.github.io/fireworks/basesite_tutorial.html?highlight=gui>`_
web gui from your local machine (assuming you have also set up
configuration files there)::

    lpad webgui



Submit the workflow
------------------------------
To launch this job through queue, use the qlaunch command from FireWorks.
qlaunch has 3 modes: singleshot, rapidfire, and multi:

* ``singleshot``: launches one job at a time
* ``rapidfire``: launches multiple jobs at once; you'll most likely
  want to use this mode where it is important to add the ``-m``
  flag to specify how many jobs to launch at once to prevent submitting
  too many jobs at once.
* ``multi``: creates one job with multiple fireworks runs

Here is an example command for launching one job from the terminal in the
same working directory as before::

    qlaunch singleshot

If you are not running your jobs through a queue, replace the
``qlaunch`` command with ``rlaunch``.

Monitor the workflow
------------------------------
If all went well, you can determine the status of your running jobs by
using the following command in the terminal::

    lpad get_fws -s RUNNING

or::

    lpad get_fws -s COMPLETED

If your job has failed, your can rerun it using the following command
(replacing ``fw_id`` with 1, which is the id of your firework, since
you only have one firework in your launchpad at this point)::

    lpad rerun_fws -i <fw_id>

Query the database for the results
--------------------------------------
Once this workflow is completed, you will see the generated Gaussian
input and output files as well as a ``run.json`` file that contains a
summary of the job in the same working directory.

Additionally, you can query the database for the results of your jobs
by using the InChI representation of the monoglyme molecule as a query
criteria:

.. code-block:: python

    from mispr.gaussian.utilities.db_utilities import get_db

    db = get_db()
    db.retrieve_run(inchi="InChI=1S/C4H10O2/c1-5-3-4-6-2/h3-4H2,1-2H3")[0]

This will return a dictionary of the results as they are saved in the
database. Alternatively, you can the view the results using MongoDB
Compass, and the generated documents from the run will appear like the
following in the ``runs`` collection of the ``gaussian`` database:

.. figure:: ../_static/document.png