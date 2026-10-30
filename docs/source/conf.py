# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Path setup
import os
import sys
sys.path.insert(0, os.path.abspath('../..')) # root of the repository, with functions/mandyocIO.py

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'post_processing'
author = 'João Pedro Macedo Silva'
release = 'v0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",           # Allows Markdown parsing
    "sphinx.ext.autodoc",    # Pulls docstrings from Python modules
    "sphinx.ext.napoleon",   # Understands the NumPy docstring format
    "sphinx.ext.autosummary", # Creates the API tables
    "sphinx.ext.intersphinx", # Link objects and classes to external docs
    'sphinx.ext.mathjax',
    'sphinx_rtd_theme',
]

# Parameters
napoleon_numpy_docstring = True
napoleon_use_param = False       # Creates a cleaner, table-like look for parameters
napoleon_use_ivar = True

autosummary_generate = True      # Creates the pages of the functions in generated/

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

myst_enable_extensions = ["dollarmath", "amsmath"]
myst_heading_anchors = 3 # anchors for the headings, used by the links of the Markdown files

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = []

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
