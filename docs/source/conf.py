# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys
import os
sys.path.insert(0, os.path.abspath("../.."))
print('sys.path: ',sys.path)
from rvdata import __version__

project = "RVdata"
copyright = "2024, BJ Fulton"
author = "BJ Fulton"
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.inheritance_diagram",
    "nbsphinx",
]

# nbsphinx configuration
nbsphinx_execute = "never"  # Don't execute notebooks during build (they have pre-computed outputs)
nbsphinx_allow_errors = True  # Continue build even if notebooks have errors
nbsphinx_prolog = """
.. raw:: html

    <style>
        .nbinput .prompt, .nboutput .prompt { display: none; }
    </style>
"""

# Prevent notebook headers from appearing in sidebar navigation
# This keeps notebooks as single entries without expanding their internal structure
toc_object_entries = False

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

missing_value_placeholder = "TBD"
