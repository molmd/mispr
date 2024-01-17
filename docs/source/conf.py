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
release = "0.0.4"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinxcontrib.mermaid",
    "sphinx.ext.coverage",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx_search.extension",
    "sphinx_tabs.tabs",
    "sphinx_togglebutton",
    "sphinx_favicon",
    "sphinx_immaterial",
    "sphinx_immaterial.apidoc.python.apigen",
]

# python_apigen_modules = {
#     "my_module": "api",
# }

mermaid_theme = {
    "theme": "dark",
    "themeVariables": {
        "primaryColor": "#BB2528",
        "primaryTextColor": "#fff",
        "primaryBorderColor": "#7C0000",
        "lineColor": "#F8B229",
        "secondaryColor": "#006100",
        "tertiaryColor": "#fff",
    },
}

autodoc_mock_imports = ["custodian", "tleap", "pymatgen", "openbabel"]

autosectionlabel_prefix_document = True
templates_path = ["_templates"]
exclude_patterns = []
sphinx_tabs_valid_builders = ["linkcheck"]
sphinx_tabs_disable_tab_closing = True
sphinx_immaterial_override_generic_admonitions = True


html_theme = "sphinx_immaterial"
html_static_path = ["_static"]
html_css_files = ["style.css"]

html_logo = "_static/logo.png"
html_title = "Materials informatics for structure-property relationships"

html_theme_options = {
    "features": [
        "content.code.annotate",
        "navigation.expand",
        "navigation.sections",
    ],  # (1)
    "palette": [
        {
            "media": "(prefers-color-scheme: light)",
            "scheme": "default",
            "primary": "red",
            "accent": "light blue",
            "toggle": {
                "icon": "material/toggle-switch-off-outline",
                "name": "Switch to dark mode",
            },
        },
        {
            "media": "(prefers-color-scheme: dark)",
            "scheme": "slate",
            "primary": "red",
            "accent": "light blue",
            "toggle": {
                "icon": "material/toggle-switch",
                "name": "Switch to light mode",
            },
        },
    ],
    "font": {
        "text": "Roboto",  # used for all the pages' text
        "code": "Roboto Mono",  # used for literal code blocks
    },
    "social": [
        {
            "icon": "fontawesome/brands/github",
            "link": "https://github.com/molmd/mispr",
            "name": "Source on github.com",
        },
        {
            "icon": "fontawesome/brands/twitter",
            "link": "https://twitter.com/molmd_group",
        },
    ],
}

source_suffix = ".rst"

master_doc = "index"

napoleon_google_docstring = True
napoleon_numpy_docstring = True
autoclass_content = "both"


def skip_params(app, what, name, obj, would_skip, options):
    if name in ("required_params", "optional_params"):
        return True
    return would_skip


def setup(app):
    app.connect("autodoc-skip-member", skip_params)
