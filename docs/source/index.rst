:orphan:

.. title:: MISPR documentation

.. module:: mispr

.. toctree::
   :hidden:

   overview
   keywords

.. toctree::
   :caption: Installation 🔧
   :hidden:
   :titlesonly:

   Overview <installation/index>
   Prerequisites <installation/dependencies>
   Configuration Files <installation/configuration>
   Running a Test Workflow <installation/test>

.. toctree::
   :caption: Workflows 🔀
   :hidden:

   workflows/basics
   workflows/supported
   workflows/tutorials
   workflows/custom

.. toctree::
   :caption: Resources 🖇️
   :hidden:

   resources/faq
   resources/resources

.. toctree::
   :caption: Code Documentation 📚
   :hidden:
   :titlesonly:

   mispr <mispr>

.. toctree::
   :caption: Development 💻
   :hidden:

   changelog
   citing
   license

##################################
MISPR |release| Documentation
##################################

MISPR is a Python library for computational materials science and contains
preset workflows for running complex hierarchical density functional
theory (DFT) and classical molecular dynamics (MD) simulations to compute
properties of materials.

.. figure:: _static/summary.jpeg
   :scale: 70%


*************
Installation
*************

.. grid:: 1 1 2 2

    .. grid-item::

        Install using `pip <https://pypi.org/project/mispr/>`__:

        .. code-block:: bash

            pip install mispr

.. important::
   Before you can start using MISPR, there are additional steps you need to follow.
   Please refer to the :doc:`installation guide <installation/index>` for complete setup instructions,
   including any dependencies or configuration files required.


*******************
Learning Resources
*******************

.. grid:: 1 1 2 2

    .. grid-item-card::
        :padding: 2

        About MISPR
        ^^^

        - :doc:`Overview  <overview>`
        - :doc:`Dependencies and prerequisites <installation/dependencies>`

    .. grid-item-card::
        :padding: 2

        Workflows
        ^^^

        - :doc:`Workflow basics  <workflows/basics>`
        - :doc:`Supported workflows <workflows/supported>`
        - :doc:`External learning resources <resources/resources>`

    .. grid-item-card::
        :padding: 2

        How-tos
        ^^^

        - :doc:`Tutorials <workflows/tutorials>`
        - :doc:`MISPR FAQ <resources/faq>`

    .. grid-item-card::
        :padding: 2

        Code documentation
        ^^^
        - :doc:`Subpackages  <mispr>`

************************************
Contributing / Reporting / Support
************************************
Contirbuting to MISPR can be in the form of:

* Requesting or adding new workflows and features
* Reporting or fixing bugs and issues
* Contributing to the documentation and/or examples

If you want to add or change something in the code, you can do this by
forking `MISPR on GitHub <https://github.com/molmd/mispr>`_ and
submitting a pull request.

If you submit a bug report, we will review it and move it to GitHub issues,
where its progress can be tracked.

For other inquiries, please contact us at rasha.atwi@stonybrook.edu.





