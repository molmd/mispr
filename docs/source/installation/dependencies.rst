===============================
Prerequisites
===============================

Conda environment
------------------------------
MISPR depends on a number of third party Python packages, and usually on
specific versions of those packages. In order not to interfere with third
party packages needed by other software on your machine or cluster, we
strongly recommend isolating MISPR in a conda environment. In the
following, we describe how to create a conda environment using
the `conda <https://docs.conda.io/projects/conda/en/latest/>`_ tool.

.. important::

   MISPR requires Python version 3.10 or higher. We have extensively tested MISPR with Python 3.10.

What you actually need to install
============================================
Not every dependency applies to every user -- it depends on which DFT
backend you plan to use. The table below summarizes what each backend
needs; the sections that follow give the exact commands.

.. list-table::
   :widths: 22 26 26
   :header-rows: 1

   * - Component
     - Gaussian user
     - **ORCA user**
   * - ``openbabel`` (conda-forge)
     - required
     - required
   * - MISPR itself (pip, from this repository)
     - required
     - required
   * - MolMD fork of ``pymatgen`` (pip, from GitHub)
     - required
     - required
   * - External program + license
     - Gaussian license, ``g16`` module
     - ORCA binaries (free academic download) + OpenMPI for parallel runs
   * - ``config.ini`` entry
     - ``gcmd`` / ``formchkcmd``
     - none (uses ``ORCA_CMD`` env variable instead)

In other words, if you only use the **ORCA backend**, your Python
environment needs just: conda-forge ``openbabel``, then ``pip install``
of MISPR (which pulls in FireWorks, custodian, pymongo, etc.
automatically) and the MolMD ``pymatgen`` fork -- plus the ORCA program
itself, which lives *outside* the conda environment (see
:ref:`installation/dependencies:ORCA backend (alternative to Gaussian)`).

Creating the environment: the exact recipe
============================================
Before creating a conda environment, ensure that Anaconda or Miniconda is installed on your system.
Most HPC clusters provide Anaconda as a loadable module. If you need to install it yourself, you can
install Miniconda by following the `official installation guide <https://docs.conda.io/projects/miniconda/en/latest/>`_.

Some MISPR dependencies must come from **conda** (conda-forge) and others
from **pip**, and the **order matters** -- conda-forge packages first, pip
packages second; see the warning below for why. The recipe is four steps,
run in exactly this order.

**Step 1 -- create and activate the environment**::

    conda create -n mispr_env python=3.10
    conda activate mispr_env

This makes an isolated environment named ``mispr_env`` with its own Python,
so nothing installed below can interfere with other software on the machine.
After activation your prompt gains a ``(mispr_env)`` prefix; every following
step assumes the environment is active. (To leave it later:
``conda deactivate``.)

**Step 2 -- compiled packages, from conda-forge FIRST**

This is the step that must come before any ``pip install``. ORCA and
Gaussian users need OpenBabel plus an explicitly pinned numpy::

    conda install -c conda-forge openbabel=3.1.1 "numpy<2"

These packages ship compiled binaries, which is both why they must come from
conda-forge (no reliable pip wheels) and why they must be installed *first*
-- they pull in conda-forge's numpy/scipy, and pip must find those already
present (see the warning below).

The ``"numpy<2"`` pin is not optional. MISPR's own dependency declaration
only says ``numpy >= 1.21.1``, which permits numpy 2.x -- but the MolMD
pymatgen fork installed in step 4 predates numpy 2 and breaks with it.
Pinning numpy at this step means every later ``pip install`` finds a
compatible numpy already present and leaves it alone.

**Step 3 -- MISPR itself, via pip**::

    pip install git+https://github.com/RuiqiLuo/mispr-psi4.git

pip sees numpy/scipy/pandas already provided by conda-forge and reuses them
instead of replacing them -- this is the payoff of the ordering. MISPR's
pure-Python dependencies (FireWorks, custodian, pymongo, mdproptools, ...)
are pulled in automatically here; you do not install them separately.

Note the source: the PyPI ``mispr`` package does **not** include the
ORCA backend, so the install must come from this repository. (For
development mode, ``git clone`` the repository first and run
``pip install -e /path/to/your/clone`` instead -- see the
:ref:`development-mode installation steps <codes-develop-mode>`.)

**Step 4 -- the MolMD fork of pymatgen**::

    pip install pymatgen@git+https://github.com/molmd/pymatgen@molmd_fix_3-9#egg=pymatgen

Required for every backend: MISPR relies on changes to pymatgen that have
not been merged upstream, so the stock pymatgen must be replaced with this
fork (see :ref:`py-package-deps` for details).

Who decides which versions get installed?
============================================
You may notice the recipe pins almost no version numbers. That is mostly
deliberate -- the versions are decided at three levels:

* **Pinned by the recipe itself** (the only things *you* control):
  Python 3.10, ``openbabel=3.1.1``, ``numpy<2``, and the exact git
  branches/tags of MISPR and the pymatgen fork.
* **Pinned by MISPR's own dependency declaration** -- when pip installs
  MISPR in step 3, it enforces MISPR's internal constraints automatically,
  most notably ``custodian==2024.10.16`` and ``pymongo<=3.12.0``. You never
  type these, and you should not "upgrade" them by hand: newer pymongo (4.x)
  in particular breaks the FireWorks version MISPR uses.
* **Left free** -- everything else (scipy, pandas, matplotlib, ...), where
  any recent release works.

The one constraint the automatic resolution does **not** protect you from is
numpy 2.x -- hence the explicit ``"numpy<2"`` pin in step 2 (see the
explanation there).

For reference, a full version set this documentation was validated against
-- useful when debugging a mysterious environment, not as targets to install
manually:

.. list-table::
   :widths: 30 30
   :header-rows: 1

   * - Package
     - Validated version
   * - Python
     - 3.10
   * - numpy
     - 1.23.5 (any 1.2x works; **not** 2.x)
   * - scipy
     - 1.15.2
   * - pandas
     - 2.3.3
   * - pymatgen (MolMD fork)
     - 2023.9.25
   * - FireWorks
     - 2.0.3
   * - custodian
     - 2024.10.16 (pinned by MISPR)
   * - pymongo
     - 3.12.0 (pinned by MISPR)
   * - openbabel
     - 3.1.1

To compare your environment against this table, run ``conda list`` inside
the activated environment (it shows conda- and pip-installed packages
together, with the installing channel in the last column).

.. warning::
   **Why the order matters -- the numpy/scipy/pandas conflict.** pip wheels
   and conda-forge builds of ``numpy``, ``scipy``, and ``pandas`` are not
   always binary-compatible with each other. If pip installs its own numpy
   first (as a dependency of mispr) and conda-forge's openbabel is
   added on top later, imports start failing with errors like::

       numpy.dtype size changed, may indicate binary incompatibility
       ValueError: All ufuncs must have type numpy.ufunc

   Installing the conda-forge packages *first* (step 2) and pip packages
   *second* (steps 3--5) avoids this entirely. Two rules keep the
   environment healthy afterwards:

   * never run ``pip install numpy`` / ``scipy`` / ``pandas`` (or
     ``pip install --upgrade`` anything that pulls them in) in this
     environment;
   * if the environment does break with the errors above, repair it with
     ``conda install -c conda-forge numpy scipy pandas --force-reinstall``
     -- or, often faster in practice, delete the environment and redo the
     recipe from step 1 in the correct order.

   **How to check who owns a package.** ``conda list numpy`` shows the
   installing channel in the last column: ``conda-forge`` means conda owns
   it (good), while ``pypi`` means pip has replaced it -- the state that
   causes the errors above if openbabel was installed from
   conda-forge. This is the first thing to check when imports start
   failing after an installation step.

.. note::
   You may need to install ``pip`` and ``setuptools`` in your conda
   environment in case the system or user version of these tools is old::

    conda install pip setuptools

Verifying the environment
============================================
After finishing the recipe, run this quick check inside the activated
environment -- each line should print without an error::

    python -c "import mispr; print('mispr', mispr.__version__)"
    python -c "from openbabel import pybel; print('openbabel OK')"
    python -c "import pymatgen; print('pymatgen', pymatgen.__version__)"
    python -c "import fireworks; print('fireworks', fireworks.__version__)"

ORCA users additionally check that the ORCA binary itself runs (it is not a
Python package -- see the ORCA section below). Executing it with no
arguments should print a usage message such as
"This program requires the name of a parameterfile as argument"::

    $ORCA_CMD    # or: /full/path/to/orca

If it instead fails with an ``error while loading shared libraries``
message, see the ORCA section below.

If any import fails with a numpy/binary-incompatibility message, see the
warning above ("How to check who owns a package").

Required Software
---------------------------------

.. _chem-software-deps:

Computational chemistry software dependencies
=============================================
At the backend, MISPR uses:

.. list-table:: 
   :widths: 20 25 40 15
   :header-rows: 1

   * - Software
     - License Type
     - Purpose
     - Installation
   * - `Gaussian <https://gaussian.com>`_
     - Commercial
     - Perform DFT calculations
     - License required
   * - `ORCA <https://www.faccts.de/orca/>`_
     - Free for academic use (registration required)
     - Perform DFT calculations (alternative to Gaussian)
     - Download from the `ORCA forum <https://orcaforum.kofo.mpg.de>`_
   * - `AmberTools24 <https://ambermd.org/AmberTools.php>`_
     - Open Source
     - Generate GAFF parameters
     - Direct download
   * - `LAMMPS <https://www.lammps.org>`_
     - Open Source
     - Run MD simulations
     - Direct download
   * - `Packmol <https://m3g.github.io/packmol/download.shtml>`_
     - Open Source
     - Build initial MD configurations
     - Follow `user guide <https://m3g.github.io/packmol/userguide.shtml>`_
   * - `Schrödinger <https://www.schrodinger.com/>`_
     - Commercial (Academic license available)
     - Generate OPLS2005 parameters (Optional)
     - License required

Ensure that you have access to the executables of these software
before using MISPR. If Gaussian, AmberTools, Schrödinger and LAMMPS are already installed on HPC
machines, the user typically needs to load their corresponding modules
before their use.

ORCA backend (alternative to Gaussian)
=======================================
As of version 0.0.5, MISPR also supports running DFT calculations through
`ORCA <https://www.faccts.de/orca/>`_ instead of Gaussian. Like Gaussian,
ORCA runs as an external program: MISPR writes a text input file,
invokes the ``orca`` binary, and parses the text output. This has two
practical consequences:

* ORCA is **not** installed into your conda environment -- it is a set of
  pre-compiled binaries you download and extract anywhere on the machine
  (no compilation needed);
* MISPR only needs to know *where* the ``orca`` executable is (see
  "Telling MISPR where ORCA is" below).

**Step 1 -- Download.** ORCA is free for academic use but requires
registration: create an account on the
`ORCA forum <https://orcaforum.kofo.mpg.de>`_, then download the Linux
archive matching your cluster's architecture (e.g.
``orca_6_1_1_linux_x86-64_shared_openmpi418.tar.xz``). The
``shared_openmpi418`` part of the name means: dynamically linked binaries
(``shared``) that expect OpenMPI 4.1.x on the machine for parallel runs.

**Step 2 -- Extract.** Transfer the archive to your machine and extract it::

    tar -xf orca_6_1_1_linux_x86-64_shared_openmpi418.tar.xz

.. warning::
   The extracted installation is large (~17 GB for ORCA 6.1). On HPC clusters
   with small home-directory quotas, extract it into a scratch/project
   filesystem instead of ``$HOME``.

.. important::
   Keep the extracted directory **intact**. The ``orca`` executable is a
   driver that calls dozens of helper binaries (``orca_scf``, ``orca_mdci``,
   ...) and shared libraries sitting next to it in the same directory. Do
   not copy the ``orca`` binary elsewhere on its own -- point MISPR at it
   inside the extracted directory.

**Step 3 -- Shared libraries.** With the ``shared`` builds, the dynamic
linker must be able to find the ``.so`` files shipped in the ORCA
directory. If running ``orca`` fails with
``error while loading shared libraries: libopenblas...``, add the
directory to ``LD_LIBRARY_PATH``::

    export LD_LIBRARY_PATH=/path/to/orca_6_1_1_linux_x86-64_shared_openmpi418:$LD_LIBRARY_PATH

**Step 4 -- OpenMPI (parallel runs only).** ORCA does not bundle MPI. To
run with ``num_cores`` > 1 you need an OpenMPI installation whose major
version matches the archive name (``openmpi418`` -> OpenMPI 4.1.x), with
``mpirun`` on your ``PATH``. On HPC clusters this is usually
``module load openmpi/4.1.x``; check with ``mpirun --version``. Serial
runs (``num_cores`` = 1, the MISPR default) work without any MPI.

Telling MISPR where ORCA is
++++++++++++++++++++++++++++++++++++++++

MISPR locates the ORCA executable through (in priority order):

1. the ``orca_cmd`` argument accepted by every ORCA workflow function;
2. the ``ORCA_CMD`` environment variable;
3. a bare ``orca`` looked up on your ``PATH``.

Note that ORCA does **not** use ``config.ini`` (that file only holds the
Gaussian/LAMMPS/AmberTools commands). The simplest setup is to export the
full path once in your ``.bashrc``/job script::

    export ORCA_CMD=/path/to/orca_6_1_1_linux_x86-64_shared_openmpi418/orca

.. note::
   Running ORCA in parallel (``num_cores`` > 1) requires ``ORCA_CMD`` (or
   ``orca_cmd``) to be the binary's **full absolute path** -- an ORCA/OpenMPI
   requirement, not a MISPR one: when ORCA re-launches itself through
   ``mpirun``, a bare ``orca`` would resolve incorrectly.

.. _py-package-deps:

Python package dependencies
=================================
* `pymatgen <https://pymatgen.org>`_: MISPR uses pymatgen for handling
  different molecule representations and i/o operations specific to
  Gaussian and LAMMPS. We have made changes to the pymatgen library to
  make it compatible with our needs in MISPR. These changes have not
  been merged yet with the main pymatgen library. Therefore, in order
  to use MISPR, you need to install the MolMD version of pymatgen
  (this is step 4 of the environment recipe above)::

    pip install pymatgen@git+https://github.com/molmd/pymatgen@molmd_fix_3-9#egg=pymatgen

* `FireWorks <https://materialsproject.github.io/fireworks/>`_: MISPR
  uses FireWorks to design, manage, and execute workflows.

  Further details can be found in the `FireWorks documentation  <https://materialsproject.github.io/fireworks/installation.html>`_.

  .. note::
   While FireWorks is used in MISPR for managing the DFT and MD
   workflows due to its many advantages, it takes some time to learn
   and get used to it.

* `custodian <https://materialsproject.github.io/custodian/>`_: MISPR uses
  custodian for handling errors that occur during the simulations and
  correcting them according to predefined rules. We have contributed a Gaussian
  plug-in to the custodian library, and these changes have been merged with 
  the main custodian library.

* `OpenBabel <https://openbabel.org>`_ to handle molecule operations 
  via pymatgen as an interface. You can install OpenBabel using conda::

    conda install -c conda-forge openbabel=3.1.1

* `MDPropTools <https://github.com/molmd/mdproptools>`_: MISPR uses mdproptools, which is a standalone 
  Python package we developed for analyzing molecular dynamics trajectories and 
  output files. 

.. note::
   FireWorks, custodian, and MDPropTools will be automatically installed as dependencies when you 
   install MISPR. You don't need to install them separately.

MongoDB
-------------------------
MISPR uses `MongoDB <https://docs.mongodb.com/manual/>`__ as the backend database.
MongoDB is a NoSQL database that is designed to store and retrieve
data in a highly efficient and scalable manner. It stores data in the
form of documents represented in the JSON (JavaScript Object Notation)
format, which is similar to a Python dictionary.

MISPR uses MongoDB to:

* Add, remove, and search the status of workflows - feature of
  `FireWorks <https://materialsproject.github.io/fireworks/>`__  (required)
* Create computational databases of DFT and MD predicted properties -
  Feature of MISPR (optional but strongly recommended)

Setting up MongoDB
============================
Options for getting MongoDB are:

* Install it yourself locally by following the instructions at
  `MongoDB <https://www.mongodb.com/docs/manual/installation/>`__.
  This is pretty simple and typically works well if you are starting out
  with MISPR and want to learn how to use a database. However, with this
  option, you are limited with the storage space on your local machine and
  you do not have the option to share the database with other users. You
  also need to have the necessary privileges to install mongo on your machine.
* Set up an account using a commercial service, which is typically
  the simplest and easiest to use but is not free of charge for databases
  with large size. Examples of such services include Atlas and MongoDB Atlas,
  which offer 500 MB databases for free. This is typically enough to get
  started for small projects.
* Self-host a MongoDB server or ask your supercomputing center to offer
  MongoDB hosting. This is more complicated than the other options and
  will require continuous maintenance of the server.

After creating a new database, you need to keep record of your credentials.
These will be used later in setting up the configuration files required
by FireWorks.

.. note::
   MongoDB must be accessible from the computers you are using to run
   the workflows.

Setting up MongoDB Atlas for an HPC cluster, step by step
==========================================================
The most common real-world setup -- a free Atlas cluster, accessed from an
HPC cluster -- involves a firewall step that trips up nearly everyone, so
here is the complete path.

**Step 1 -- create the cluster.** Sign up at
`MongoDB Atlas <https://www.mongodb.com/cloud/atlas>`_ and create a free
(M0, 512 MB) cluster. The Atlas *account* is for the web dashboard only; it
is not what your scripts will log in with. Note that the setup wizard
silently whitelists the IP address *you are browsing from* ("Add My Current
IP Address") -- which is why connections sometimes "just work" from
machines on the same network as your browser, while failing from the HPC
cluster. Step 3 below is about fixing exactly that.

**Step 2 -- create a database user.** In the dashboard under
*Security -> Database Access*, add a database user with a username and
password (role: "Read and write to any database" is fine to start). This
user -- not your Atlas account -- is what goes into the connection string.
Prefer a password without special characters (letters/digits only): special
characters must be URL-encoded inside a connection string (``@`` becomes
``%40``, etc.), a classic source of "Authentication failed" confusion.

**Step 3 -- open the firewall (Network Access).** This is the step where
most people get stuck. By default **Atlas rejects every connection**, from
everywhere; under *Security -> Network Access* you must whitelist the IP
addresses your connections will come *from* -- i.e. the outward-facing IP
of the HPC cluster, not of your laptop (and not the cluster's internal
``10.x``/``192.168.x`` addresses). Find it by running, on a login node::

    curl -s ifconfig.me

and add that IP in *Network Access -> Add IP Address*. Two HPC realities to
plan around:

* **Compute nodes may leave through different IPs** than the login node (or
  through a whole NAT range). If a connection test works on the login node
  but workflows fail inside jobs, this is why. Options, in decreasing order
  of niceness: ask your HPC support for the cluster's outbound IP range and
  whitelist that; or whitelist ``0.0.0.0/0`` ("allow access from
  anywhere"), which sounds alarming but is common practice for student/lab
  Atlas clusters -- the connection is still TLS-encrypted and
  password-protected, the firewall is simply no longer an *additional*
  layer. Check your institution's policy if in doubt.
* **Some clusters' compute nodes have no internet access at all** (only
  login nodes do). No Atlas whitelist can fix that -- test early (step 5)
  and, if so, ask your HPC support whether an outbound proxy exists or
  whether MongoDB can be hosted inside the cluster instead.

**Step 4 -- get the connection string.** Dashboard -> *Connect* ->
*Drivers* -> Python. You get an SRV-style URI::

    mongodb+srv://<user>:<password>@cluster0.xxxxx.mongodb.net/?retryWrites=true&w=majority

Insert the name of the database MISPR should write to (any name you like,
e.g. ``mispr``) as a path before the ``?``::

    mongodb+srv://<user>:<password>@cluster0.xxxxx.mongodb.net/mispr?retryWrites=true&w=majority

MISPR and FireWorks both parse the database name out of that path, so the
URI is the *only* place it needs to be stated. (SRV URIs need the
``dnspython`` package -- installed automatically with MISPR.)

**Step 5 -- test from the cluster.** In your activated environment **on the
HPC cluster** (not your laptop -- the whitelist is per-IP!)::

    python -c "
    from pymongo import MongoClient
    client = MongoClient('mongodb+srv://<user>:<password>@cluster0.xxxxx.mongodb.net/mispr?retryWrites=true&w=majority')
    print(client.list_database_names())
    "

A list of database names printed = everything works. The two common
failures map directly to the steps above:

* ``ServerSelectionTimeoutError`` mentioning ``SSL handshake failed ...
  tlsv1 alert internal error`` -- the *signature* symptom of a
  non-whitelisted IP: Atlas actively cuts the connection during the TLS
  handshake. Whitelist the IP of the node you are on (step 3) -- and note
  that different login/compute nodes of the same cluster often leave
  through *different* IPs, so prefer whitelisting the whole range (e.g.
  ``129.49.82.0/24``) over a single address.
* ``ServerSelectionTimeoutError`` with plain timeouts (no SSL mention) --
  the connection never reached Atlas at all: the node has no internet
  access, or an outbound firewall blocks port 27017.
* ``Authentication failed`` -- the connection reached Atlas but the
  credentials were rejected: wrong database user (step 2), or special
  characters in the password that were not URL-encoded.

Then repeat the same test inside a short compute-node job (``srun``/an
``sbatch`` script) -- this catches the compute-nodes-different-IP problem
of step 3 *before* your first real workflow mysteriously hangs.

**Step 6 -- put the URI into the configuration files.** Both files accept
the URI directly; see :doc:`Configuration Files <configuration>` for the
full file layouts. In ``db.json``, set ``uri_mode`` to true and put the
full URI in ``host`` (the other credential fields are then not used)::

    {
        "host": "mongodb+srv://<user>:<password>@cluster0.xxxxx.mongodb.net/mispr?retryWrites=true&w=majority",
        "uri_mode": true,
        "admin_user": null,
        "admin_password": null,
        "database": null,
        "collection": null,
        "aliases": {}
    }

In ``my_launchpad.yaml``, set ``uri_mode: true`` and put a URI in ``host``
the same way -- but give it a *different* database name in the path (e.g.
``.../fireworks?retryWrites=...``), so FireWorks' internal bookkeeping and
MISPR's results live in separate databases on the same cluster.

Testing a local (non-Atlas) MongoDB
================================================
If you instead installed MongoDB yourself on the machine you run on, the
connection test is the same idea with the default local address::

    python -c "
    from pymongo import MongoClient
    client = MongoClient('mongodb://localhost:27017/')
    print(client.list_database_names())
    "

If MongoDB is running on a different machine or port, adjust the connection
string accordingly.