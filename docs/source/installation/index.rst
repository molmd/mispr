=============
Installation
=============
Before installing MISPR, you need to follow the steps below in order:

1. (Optional) Create a :ref:`installation/dependencies:Conda environment`
2. Make sure you have access to the :ref:`computational chemistry software dependencies <chem-software-deps>`
   needed to run the DFT and MD simulations
3. Install :ref:`python package dependencies <py-package-deps>`
4. Set up :ref:`installation/dependencies:MongoDB` database
5. :ref:`Install MISPR <installation/index:Installing MISPR>`
6. Prepare the :doc:`configuration files <configuration>`
7. :doc:`Run a test workflow <test>`

.. note::
   For the DFT engine (step 2), a Gaussian license is **not** the only
   option: the ESP, bond dissociation energy, binding energy, and redox
   potential workflows can alternatively run through
   `ORCA <https://www.faccts.de/orca/>`_ (free for academic use). See the
   corresponding section in :doc:`Prerequisites <dependencies>` for how to
   install it, and :doc:`Workflow Tutorials <../workflows/tutorials>` for
   usage.

.. note::
   Throughout the installation instructions, it is assumed that you are
   familiar with Python and with basic Linux shell commands. If not,
   `Linux Journey <https://linuxjourney.com/lesson/the-shell>`_ and
   `Python For Beginners <https://www.python.org/about/gettingstarted/>`_
   are some recommended starting points.

Installing MISPR
--------------------------------
MISPR can be installed either from the python package
index (good for most users) or directly from its GitHub
repository (good for developers).

Installation Method 1: Using pip
================================
For most users, install directly from this repository:

.. code-block:: bash

    pip install git+https://github.com/RuiqiLuo/mispr-psi4.git

.. warning::
   Do not use ``pip install mispr`` (the PyPI package) with this
   documentation: the PyPI release is the upstream version and does **not**
   include the ORCA backend described here.

Installation Method 2: Development mode
=======================================

.. _codes-develop-mode:

For developers or users who need to modify the source code, install MISPR in development mode. 

.. note::
   If you had already installed MISPR via pip or conda, you
   should uninstall that first before starting the installation in
   development mode. This ensures that you will not have any conflicts
   resulting from two different code installations.

The steps for installing the package in development mode are below.

1. Activate your conda environment or virtual environment

2. Create a ``codes`` directory wherever you keep your projects
   (e.g. ``~/codes``; on HPC clusters with small home quotas, a
   scratch/project filesystem is a better choice)

3. ``cd`` into your newly created ``codes`` directory

4. Clone the package you want to install in development mode using git::

    git clone https://github.com/RuiqiLuo/mispr-psi4.git

   (This fork contains the ORCA backend documented here; the
   upstream repository is `molmd/mispr <https://github.com/molmd/mispr>`_.)
   Now you should have the repository directory in your ``codes`` directory.

5. ``cd`` into the mispr directory and run
   ``pip install -e .`` or use the ``conda`` equivalent. Once installed,
   if you make changes to the code, the changes
   will take effect immediately without having to reinstall the package.

Post-installation
-------------------------
1. Before you go any further, confirm your package installations are correct::

    ipython -c "import mispr; print('mispr', mispr.__version__)"

   This should print the version number without any errors. For a fuller
   check (openbabel, pymatgen, FireWorks, and -- for ORCA users -- the
   ORCA binary itself) plus the exact package versions this documentation
   was validated against, see
   :ref:`installation/dependencies:Verifying the environment`.

2. To update the mispr code later on: in development mode, execute
   ``git pull`` in the repository directory (with ``pip install -e .`` the
   changes take effect immediately -- but running Jupyter kernels must be
   restarted to pick them up). If you installed via pip from the repository
   URL, re-run
   ``pip install --force-reinstall --no-deps git+https://github.com/RuiqiLuo/mispr-psi4.git``.
   Do **not** use ``pip install --upgrade mispr`` -- that pulls the PyPI
   release, which would replace this version with one lacking the ORCA
   backend.