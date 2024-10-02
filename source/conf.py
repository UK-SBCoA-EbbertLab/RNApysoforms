# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

## IMPORT STATEMENTS
#import sphinx_rtd_theme  # Read the Docs theme
#import os
#import sys

# Ensure that Sphinx can find the 'src' directory where your code is located
#sys.path.insert(0, os.path.abspath('../src/'))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pytranscript'
copyright = '2024, Bernardo Aguzzoli Heberle'
author = 'Bernardo Aguzzoli Heberle'
release = '0.1.0-dev'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',          # Automatically document your code
    'sphinx.ext.napoleon',         # Support for NumPy and Google style docstrings
    'sphinx.ext.viewcode',         # Add links to highlighted source code
    'sphinx.ext.autosummary',      # Generate summary pages for modules/classes
    'sphinx.ext.autodoc.typehints' # Automatically include type hints in docs (if using type annotations)
]

# Path to templates used in the documentation
templates_path = ['_templates']

# Patterns to exclude when looking for source files
exclude_patterns = []

# Automatically generate autosummary pages
autosummary_generate = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Set the theme to Read the Docs
html_theme = 'sphinx_rtd_theme'

# Path to custom static files, such as CSS or JavaScript
html_static_path = ['_static']
