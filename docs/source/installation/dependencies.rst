===============================
Prerequisites
===============================

Virtual python environment
------------------------------
MISPR depends on a number of third party Python packages, and usually on
specific versions of those packages. In order not to interfere with third
party packages needed by other software on your machine or cluster, we
strongly recommend isolating MISPR in a virtual python environment. In the
following, we describe how to create a virtual python environment using
the `virtualenv <https://virtualenv.pypa.io/en/latest/>`_ tool, but
feel free to use your preferred environment manager
(e.g. `conda <https://conda.io/docs/>`_).

.. tip::
   :title: Recommendation

   We recommend that you use Python version 3.7 since this is the version that we have
   tested MISPR with.

Creating the virtual environment
=================================
To create and activate a new virtual environment, go to your
``|CODES_DIR|`` (see :doc:`Definition <../keywords>`), and run the following commands::

    pip install --user --upgrade virtualenv     # install virtualenv tool
    python -m virtualenv mispr_env              # create "mispr_env" environment
    source mispr_env/bin/activate               # activate "mispr_env" environment

This will create a directory in your ``|CODES_DIR|`` named ``mispr_env``,
where all the packages will be installed. After activation, your prompt
should have ``(mispr_env)`` in front of it, indicating that you are
working inside the virtual environment. The activation script ensures
that python programs have access only to packages installed inside the
virtualenv.
To deactivate the enviornment, simply run::

    deactivate

.. note::
   You may need to install ``pip`` and ``setuptools`` in your virtual
   enviornment in case the system or user version of these tools is old::

    pip install -U setuptools pip

Computational chemistry software
---------------------------------

At the backend, MISPR uses:

* `Gaussian <https://gaussian.com>`_ software to perform DFT calculations
* `AmberTools <https://ambermd.org/AmberTools.php>`_  to generate GAFF parameters
* `LAMMPS <https://www.lammps.org/#gsc.tab=0>`_ to run MD simulations
* `Packmol <https://m3g.github.io/packmol/download.shtml>`_ to
  create initial configurations for MD simulations. To install packmol,
  follow their `user guide <https://m3g.github.io/packmol/userguide.shtml>`_

Ensure that you have access to the executables of these software
before using MISPR. Gaussian is a commercial software
that requires a license while AmberTools, LAMMPS, and Packmol are open source.
If Gaussian, AmberTools, and LAMMPS are already installed on supercomputing
resources, the user typically needs to load their corresponding modules
before their use.

Materials Project base libraries
---------------------------------
* `pymatgen <https://pymatgen.org>`_: MISPR uses pymatgen for handling
  different molecule representations and i/o operations specific to
  Gaussian and LAMMPS. We have made changes to the pymatgen library to
  make it compatible with our needs in MISPR. These changes have not
  been merged yet with the main pymatgen library. Therefore, in order
  to use MISPR, you need to install the MolMD version of pymatgen by
  running the following commands in your ``|CODES_DIR|``::

    git clone https://github.com/molmd/pymatgen.git
    cd pymatgen
    python setup.py install
* `FireWorks <https://materialsproject.github.io/fireworks/>`_: MISPR
  uses FireWorks to design, manage, and execute workflows. To install,
  simply type::

    pip install FireWorks

  Further details can be found in the `FireWorks documentation  <https://materialsproject.github.io/fireworks/installation.html>`_.

  .. note::
   While FireWorks is used in MISPR for managing the DFT and MD
   workflows due to its many advantages, it takes some time to learn
   and get used to it.

* `custodian <https://materialsproject.github.io/custodian/>`_: MISPR uses
  custodian for handling errors that occur during the simulations and
  correcting them according to predefined rules. We have added a Gaussian
  plug-in to the custodian library, but similar to the pymatgen changes,
  these changes have not been merged yet with the main custodian library.
  Therefore, in order to use MISPR, you need to install the MolMD version
  of custodian by running the following commands in your ``|CODES_DIR|``::

    git clone https://github.com/molmd/custodian.git
    cd custodian
    python setup.py install

MongoDB
-------------------------
Following the design decisions of the Materials Project, MISPR uses
`MongoDB <https://docs.mongodb.com/manual/>`__ as the backend database.
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

Testing your MongoDB connection
================================
**Establishing a Connection to MongoDB Using Pymongo:**

You need to import MongoClient from pymongo and then create a new MongoClient instance.
This instance is used to connect to your MongoDB instance:

.. code-block:: python

    from pymongo import MongoClient

    client = MongoClient("mongodb://localhost:27017/")

In this example, we're connecting to a MongoDB instance that runs on the same machine
(localhost) on port 27017, which is the default port for MongoDB.

**Testing the Connection to MongoDB:**

We can check the connection by listing all the databases:

.. code-block:: python

    print(client.list_database_names())

If the connection is successful, this command will return a list of names of the databases that are present in the
MongoDB instance.

Remember, for you to connect to a MongoDB instance, the MongoDB server needs to be installed and running.
If it's not running on localhost:27017, you will need to provide the appropriate connection string.