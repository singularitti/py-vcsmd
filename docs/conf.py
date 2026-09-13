"""Sphinx configuration for the VCSMD documentation."""

from importlib.metadata import version as package_version

project = "VCSMD"
release = package_version("vcsmd")
version = release
language = "en"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinxcontrib.mermaid",
]
exclude_patterns = ["_build", ".DS_Store", "Thumbs.db"]

myst_enable_extensions = ["amsmath", "colon_fence", "dollarmath"]
myst_fence_as_directive = ["mermaid"]
myst_heading_anchors = 3
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_preserve_defaults = True

html_theme = "furo"
html_title = f"VCSMD {release}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sourcelink = False
html_copy_source = False
html_show_copyright = False
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#087f8c",
        "color-brand-content": "#086c78",
    },
    "dark_css_variables": {
        "color-brand-primary": "#6edbe5",
        "color-brand-content": "#6edbe5",
    },
}

# Pin browser-side renderers independently of the Python dependency ranges.
mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js"
mermaid_version = "11.12.1"
mermaid_init_config = {"startOnLoad": False, "securityLevel": "strict"}
