# -- Project information -----------------------------------------------------
project = 'RNApysoforms'
copyright = '2024, Bernardo Aguzzoli Heberle'
author = 'Bernardo Aguzzoli Heberle'
release = '0.1.0-dev'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx_plotly_directive',
    'nbsphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc.typehints',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = ['_build', '**.ipynb_checkpoints']
autosummary_generate = True
add_module_names = True
nbsphinx_execute = 'never'

# Plotly configuration
import plotly.io as pio
pio.renderers.default = 'sphinx_gallery'

# nbsphinx configuration
nbsphinx_requirejs_path = ''
nbsphinx_widgets_language = 'javascript'

# HTML output configuration
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_js_files = [
    'https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js',
]


# Change the max-width of the content area
html_theme_options = {
    'body_max_width': '80%',  # Adjust the width as needed
}

# Add custom CSS file
def setup(app):
    app.add_css_file('custom.css')

# Specify source suffixes
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}