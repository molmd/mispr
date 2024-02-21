=============
Installation
=============
Before installing MISPR, you need to follow the steps below in order:

1. (Optional) Create a :ref:`installation/dependencies:Virtual python environment`
2. Make sure you have access to the :ref:`installation/dependencies:Computational chemistry software`
   needed to run the DFT and MD simulations
3. Install :ref:`installation/dependencies:Materials Project base libraries`
4. Set up :ref:`installation/dependencies:MongoDB` database
5. :ref:`Install MISPR and MDPropTools <installation/index:Installing MISPR and MDPropTools>`
6. Prepare the :doc:`configuration files <configuration>`
7. :doc:`Run a test workflow <test>`

.. note::
   Throughout the installation instructions, it is assumed that you are
   familiar with Python and with basic Linux shell commands. If not,
   `Linux Journey <https://linuxjourney.com/lesson/the-shell>`_ and
   `Python For Beginners <https://www.python.org/about/gettingstarted/>`_
   are some recommended starting points.

Installing MISPR and MDPropTools
--------------------------------
MISPR and MDPropTools can be installed either from the python package
index (good for most users) or directly from their GitHub
repositories (good for developers).

Installation Method 1: Using pip
================================
To install, simply type:

.. code-block:: bash

    pip install mispr
    pip install mdproptools

Installation Method 2: Development mode
=======================================

.. _codes-develop-mode:

To make changes directly to the source and contribute to the development
of MISPR, you can install MISPR and MDPropTools in development mode.

.. note::
   If you had already installed MISPR via pip or conda, you
   should uninstall that first before starting the installation in
   development mode. This ensures that you will not have any conflicts
   resulting from two different code installations.

The steps for installing the packages in development mode are below.

1. Activate your conda environment or virtual environment

2. Create a ``codes`` directory in ``|CODES_DIR|``

3. ``cd`` to your newly created ``|CODES_DIR|/codes`` directory

4. Clone the packages you want to install in development mode using git::

    git clone https://github.com/molmd/mdproptools.git
    git clone https://github.com/molmd/mispr.git

   Now you should have mdproptools and mispr directories in your ``codes``
   directory.

5. For each of these packages, ``cd`` into their folders and run
   ``pip install -e .`` or use the ``conda`` equivalent. Once installed,
   if you make changes to the code in these packages, the changes
   will take effect immediately without having to reinstall the package.

Post-installation
-------------------------
1. Before you go any further, confirm your package installations are correct.
   First start IPython by typing ``ipython`` in your terminal, then confirm that
   the commands ``import mdproptools`` and ``import mispr`` execute
   without any errors

2. To update these codes later on, execute ``git pull`` followed by
   ``pip install -e .`` or the ``conda`` equivalent in the corresponding
   folder if you installed in development mode. If you installed via pip,
   you can simply execute ``pip install --upgrade mispr``.