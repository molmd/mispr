.. title:: MISPR developer installation
.. _developer installation:

=======
Installing MISPR in development mode
=======

Installation
============

.. _codes-develop-mode:

To make changes directly to the source and contribute to the development
of MISPR, you can install MISPR in development mode.

**Note**: If you had already installed MISPR via pip or conda, you
should uninstall that first before starting the installation in
development mode. This ensures that you will not have any conflicts
resulting from two different code installations.

The steps for installing the packages in development mode are below.

1. Activate your conda environment or virtual environment

2. Create a ``codes`` directory in ``|CODES_DIR|``

3. ``cd`` to your newly created ``|CODES_DIR|/codes`` directory

4. Clone the packages you want to install in development mode using git::

    git clone https://github.com/molmd/pymatgen/tree/molmd_fix.git
    git clone https://github.com/molmd/custodian.git
    git clone https://github.com/molmd/mdproptools.git
    git clone https://github.com/molmd/mispr.git

   Now you should have pymatgen, custodian, mdproptools, and mispr
   directories in your ``codes`` directory.

5. For each of these packages, ``cd`` into their folders and run
   ``pip install -e .`` or use the ``conda`` equivalent. Once installed,
   if you make changes to the code in these packages, the changes
   will take effect immediately without having to reinstall the package.

Post-installation
============
1. Before you go any further, confirm your package installations are correct.
   First start IPython by typing ``ipython`` in your terminal, then confirm that
   the commands ``import pymatgen``, ``import custodian``, ``import mdproptools``,
   and ``import mispr`` execute without any errors

2. To update these codes later on, execute ``git pull`` followed by
   ``pip install -e .`` or the ``conda`` equivalent in the corresponding
   folder


