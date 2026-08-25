===============================
Configuration Files
===============================
After setting up the environment and installing the software, you will
need to configure the software to work with your system. This is done by
creating the following set of configuration files.

.. note::
   This page is intended to help you get set-up for the first time using
   FireWorks and MISPR so you can learn how these software work. Please
   refer to the `FireWorks <https://materialsproject.github.io/fireworks/>`__
   documentation for more information on these files and how FireWorks works.
   Here, we will only discuss the basic configuration options which are
   sufficient for using MISPR as in this tutorial.


Writing the Configuration Files
------------------------------------

.. tab-set::

    .. tab-item:: db.json

        This file contains the basic MongoDB information like the
        credentials required to connect to the database where the
        calculation outputs will be stored. Note that
        JSON strings require double quotes except for the value of
        "port" which is an integer.

        .. code-block:: json

            {
                "admin_user": "|ADMIN_USERNAME|",
                "admin_password": "|ADMIN_PASSWORD|",
                "aliases": {},
                "collection": "|COLLECTION|",
                "database": "|DB_NAME|",
                "host": "|HOSTNAME|",
                "port": |PORT|,
            }

    .. tab-item:: my_fworker.yaml

        This file stores your FireWorker's credentials. In
        the `FireWorks <https://materialsproject.github.io/fireworks/index.html#centralized-server-and-worker-model>`__, a
        FireWorker can be as simple as the workstation used to host the
        LaunchPad or complicated like a supercomputing center with a
        queueing system.

        .. code-block:: yaml

            name: |WORKER_NAME|
            category: ''
            query: '{}'
            env:
                db_file: |CODES_DIR|/config/db.json
                scratch_dir: null

        The following parameters are defined in the file:

        * ``name``: the name of the worker where your job will be run;
          this is helpful when you have multiple workers; see
          `FireWorks documentation on controlling the Worker <https://materialsproject.github.io/fireworks/controlworker.html?highlight=category>`_
          if you need more information on setting up this file if you are
          using more than one worker.

        * ``category`` and ``queue``: these parameters can control which
          calculations are run on which worker; the default settings will
          allow all calculations to be run

        * ``env``: defines worker-specific settings like the path to the
          db file and the scratch directory for fast disk access

    .. tab-item:: my_launchpad.yaml

        This is the FireWorks LaunchPad file that contains the MongoDB
        credentials required to connect to the database for storing
        and managing workflows within FireWorks. Note that the ``db.json``
        file we created earlier is used to connect to the database
        where the results are stored and is used by MISPR while
        ``my_launchpad.yaml`` is used by FireWorks. The two databases
        can be the same or different databases. If they are the same databases,
        then the information here will be mostly the same as that in the
        ``db.json`` file.

        .. code-block:: yaml

                host: |HOSTNAME|
                port: |PORT|
                name: |LAUNCHPAD_NAME|
                username: |ADMIN_USERNAME|
                password: |ADMIN_PASSWORD|
                logdir: null
                strm_lvl: INFO
                user_indices: []
                wf_user_indices: []
                authsource: null
                uri_mode: |URI_MODE|
                mongoclient_kwargs: {}

        The following parameters need to be defined in the file:

        * ``host``: the hostname of the MongoDB server

        * ``port``: the port number of the MongoDB server

        * ``name``: the name of the MongoDB server

        * ``username``: the username to connect to the MongoDB server

        * ``password``: the password to connect to the MongoDB server

        Note that if the ``uri_mode`` is set to true, the ``host``
        should be the full `URI string <https://www.mongodb.com/docs/manual/reference/connection-string/>`_.
        In this case, the ``username`` and ``password`` are not used.

        If you want to pass other custom keyword arguments
        (e.g., SSL/TLS arguments) to the MongoClient connection, you
        can do that via ``mongoclient_kwargs``. See
        `pymongo documentation <https://pymongo.readthedocs.io/en/stable/api/pymongo/mongo_client.html>`_
        for more details.

    .. tab-item:: my_qadapter.yaml

        This is the queue adapter file required by FireWorks to
        automatically communicate with the queueing system.
        The example provided here is for SLURM machines and does not
        include a full list of possible parameters, but you can
        check the rest of the parameters or parameters that can be
        specified for other queue systems (e.g., PBS, SGE, etc.)
        `here <https://github.com/materialsproject/fireworks/tree/main/fireworks/user_objects/queue_adapters>`_.

        .. code-block:: yaml

                _fw_name: CommonAdapter
                _fw_q_type: SLURM
                rocket_launch: rlaunch -w |CODES_DIR|/config singleshot
                nodes: 1
                walltime: 24:00:00
                queue: null
                account: null
                job_name: null
                pre_rocket: null
                post_rocket: null
                logdir: |CODES_DIR|/logs

        The following parameters are defined in the file:

        * ``_fw_name``: ``CommonAdapter`` means that the queue is one of
          the built-in queue systems

        * ``_fw_q_type``: the queue system type (e.g., SLURM, PBS, SGE, etc.)

        * ``rocket_launch``: the method to use for launching Rockets

        * ``nodes``, ``walltime``, ``queue``, ``account``, ``job_name``:
          parameters you normally specify in your SLURM script for
          allocating resources

        * ``pre_rocket`` and  ``post_rocket``: the commands to run
          before and after launching the Rocket (e.g., module load
          packages)

        * ``logdir``: path to the log directory

        .. note::
            Specifying singleshot in the file will limit each
            reserved job to running only one firework at a time even if other
            fireworks are waiting to be run. This can be changed to rapidfire
            to run all fireworks in parallel. You can go over the FireWorks
            documentation to learn the difference between these launching modes.

    .. tab-item:: config.ini

        This file contains the commands to run Gaussian, LAMMPS, and AmberTools.
        These commands are specific to your computing resources you are
        running on. The example provided here is meant to show how these
        commands should be defined, but you need to change them to match your
        system.

        .. code-block:: yaml

                [RunCalc]
                gcmd: g16 < "$input_path$" > "$output_path$"
                formchkcmd: formchk "$input_path$" "$output_path$"

                [LammpsRunCalc]
                lcmd: mpirun -np $SLURM_NTASKS lmp_mpi -in $control_path$
                lammps_gpu_cmd: null

                [AmbertoolsRunCalc]
                acmd: antechamber -i $input_file$ -fi $input_type$ -o $output_file$ -fo $output_type$ -c $charge_method$ -s 2
                pcmd: parmchk2 -i $input_file$ -f mol2 -o $output_file$
                tcmd: tleap -f $input_file$

                [MaestroCalc]
                mae_cmd: $SCHRODINGER/utilities/structconvert $input_file$ $output_file$
                ffld_cmd: $SCHRODINGER/utilities/ffld_server -imae $input_file$ -version 14 -print_parameters -out_file $output_file$

        The following commands are defined in the file:

        * ``gcmd``: the command to run Gaussian
        * ``formchkcmd``: the command to run Gaussian formchk to convert
          a Gaussian checkpoint file into formatted forms
        * ``lcmd``: the command to run LAMMPS
        * ``lammps_gpu_cmd``: the command to run LAMMPS on a GPU
        * ``acmd``: the command to run Antechamber
        * ``pcmd``: the command to run Parmchk2
        * ``tcmd``: the command to run tleap

        .. note::
            Anything between two dollar signs ($ $) is a placeholder for
            a variable and should not be changed.

            Anything between the square brackets ([]), e.g., [RunCalc],
            or before the colons (:), e.g., gmcd, should not be changed
            since these are used to point MISPR to the commands to run.

    .. tab-item:: FW_config.yaml

        This is the master FireWorks configuration file that controls
        FireWorks settings and points to the location of the other
        configuration files.

        .. code-block:: yaml

                CONFIG_FILE_DIR: |CODES_DIR|/config

        The ``CONFIG_FILE_DIR`` is expected to contain the
        other configuration files. For a list of control settings that
        can be added to this file, check
        `FireWorks documentation on modifying the FW config <https://materialsproject.github.io/fireworks/config_tutorial.html>`_.

Configuring Bash Profile
------------------------------
After creating the above six configuration files and replacing the
placeholders with your specific settings, create a directory in
your ``|CODES_DIR|`` (see :doc:`Definition <../keywords>`) called ``config``
and move the above configuration files into it. The ``|CODES_DIR|/config``
should look like:

::

    config
    ├── config.ini
    ├── db.json
    ├── FW_config.yaml
    ├── my_fworker.yaml
    ├── my_launchpad.yaml
    └── my_qadapter.yaml

Now, append the following lines to your ``.bash_profile`` or ``.bashrc``
file in order to set an environment variable that tells FireWorks where
to find the ``FW_config.yaml`` file, which will in turn tell FireWorks
where the rest of the configuration files are:

.. code-block:: bash

    export FW_CONFIG_FILE=|CODES_DIR|/config/FW_config.yaml