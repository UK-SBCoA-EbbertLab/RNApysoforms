# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'RNApysoforms'
copyright = '2024, Bernardo Aguzzoli Heberle'
author = 'Bernardo Aguzzoli Heberle'
release = '0.1.0-dev'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'nbsphinx',                    # Jupyter notebooks vignettes
    'sphinx.ext.autodoc',          # Automatically document your code
    'sphinx.ext.napoleon',         # Support for NumPy and Google style docstrings
    'sphinx.ext.viewcode',         # Add links to highlighted source code
    'sphinx.ext.autosummary',      # Generate summary pages for modules/classes
    'sphinx.ext.autodoc.typehints' # Automatically include type hints in docs (if using type annotations)
]

# Path to templates used in the documentation
templates_path = ['_templates']

# Patterns to exclude when looking for source files
exclude_patterns = ['_build', '**.ipynb_checkpoints']

# Automatically generate autosummary pages
autosummary_generate = True
add_module_names = True  # Display only function names without the module prefix

# Don't re-run notebooks when building
nbsphinx_execute = 'never'

# Add Markdown support
extensions.append('myst_parser')


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Set the theme to Read the Docs
html_theme = 'sphinx_rtd_theme'

# Path to custom static files, such as CSS or JavaScript
html_static_path = ['_static']

## Allow plotly figures to render
def setup(app):
    app.add_js_file('https://cdn.plot.ly/plotly-latest.min.js')
