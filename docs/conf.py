# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import datetime
#sys.path.insert(0, os.path.abspath('.'))

sys.path.insert(0, os.path.abspath('../'))


# -- Project information -----------------------------------------------------
year = datetime.datetime.now().date().year
version_file = os.path.join(os.path.dirname(__file__), os.pardir, "actin2", "VERSION")

try:
    with open(version_file, 'r') as file:
        version = file.read()
except FileNotFoundError:
    version = "unknown"

project = 'ACTIN'
copyright = f'2018-{year}, J. Gomes da Silva'
author = 'J. Gomes da Silva'

# The full version, including alpha/beta/rc tags
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
#import sphinx_rtd_theme
#import myst_parser
#import nbsphinx

extensions = [
    #"myst_parser",
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    'sphinx.ext.napoleon',
    'nbsphinx',
]

napoleon_google_docstring = True
napoleon_numpy_docstring = True

napoleon_include_special_with_doc = True
napoleon_include_private_with_doc = True

#pygments_style = 'sphinx'
html_show_sourcelink = True

html_theme = "sphinx_rtd_theme"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
#html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = 'img/logo.png'
