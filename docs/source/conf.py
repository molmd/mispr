# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.dirname(".."))
sys.path.insert(0, os.path.dirname("../../mispr"))
sys.path.insert(0, os.path.abspath("../.."))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "MISPR"
copyright = "2022, MolMD Group"
author = "MolMD Group"
release = "0.0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosectionlabel",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ["_static"]

html_theme_options = {
    "canonical_url": "https://molmd.github.io/mispr",
    "logo_only": True,
    "display_version": True,
    "prev_next_buttons_location": None,
    "style_external_links": True,
    "style_nav_header_background": "linear-gradient(0deg, rgba(23,63,162,1) 0%, rgba(0,70,192,1) 100%)",
    "collapse_navigation": True,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}

html_context = {
    "display_github": True,
    "github_user": "molmd",
    "github_repo": "mispr",
    "github_version": "master",
    "conf_py_path": "/docs/source/",
}

master_doc = "index"