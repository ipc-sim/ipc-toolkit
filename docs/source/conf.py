# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import subprocess
import pathlib

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

import ipctk

project = "IPC Toolkit"
copyright = '2020-2023, IPC-Sim Organization; MIT License'
author = "Zachary Ferguson"
version = ipctk.__version__

# -- General configuration ---------------------------------------------------

# Doxygen
pathlib.Path("../build/doxyoutput").mkdir(parents=True, exist_ok=True)
if (not subprocess.run(["doxygen", "Doxyfile"])):
    print("Doxygen failed! Exiting")
    exit(1)

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "autoclasstoc",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    # "sphinx.ext.autosectionlabel",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.inheritance_diagram",
    "sphinx.ext.graphviz",
    "breathe",
    "myst_parser",
    "nbsphinx",
    "sphinx_immaterial",
    "sphinx_immaterial.apidoc.python.apigen",
    "sphinx_immaterial.apidoc.format_signatures",
    # 'sphinx_autodoc_toolbox.collapse',
    "sphinxcontrib.bibtex",
    "sphinxemoji.sphinxemoji",
    "sphinx_last_updated_by_git",
]

bibtex_bibfiles = ['refs.bib']
bibtex_reference_style = 'author_year'
bibtex_default_style = 'plain'

myst_enable_extensions = [
    "dollarmath",
]

object_description_options = [
    ("cpp:.*", dict(clang_format_style={"BasedOnStyle": "WebKit"})),
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

root_doc = "index"

suppress_warnings = ["myst.header"]

# Setup the breathe extension
breathe_projects = {
    project: "../build/doxyoutput/xml"
}
breathe_default_project = project
breathe_default_members = (
    "members",
    "undoc-members",
    "protected-members",
    "private-members",
)
breathe_show_define_initializer = True
# breathe_show_include = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": True,
    'special-members': True,
    'show-inheritance': True,
}

# -- GraphViz configuration ----------------------------------
graphviz_output_format = 'svg'

graphviz_dot_args = ["-Ecolor=#CE93D8", "-Kdot"]

# python_apigen_modules = {
#     "ipctk": "",
# }

# python_apigen_default_groups = [
#     ("class:.*", "Classes"),
#     (r".*\.__(init|new)__", "Constructors"),
#     (r".*\.__(str|repr)__", "String representation"),
# ]

# Tell sphinx what the primary language being documented is.
primary_domain = "cpp"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "cpp"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "sphinx_immaterial"

# Material theme options
html_theme_options = {
    "palette": [{
        "media": "(prefers-color-scheme: light)",
        "scheme": "default",
        "primary": "deep-purple",
        "accent": "deep-purple",
        "toggle": {
            "icon": "material/brightness-7",
            "name": "Switch to dark mode",
        }
    }, {
        "media": "(prefers-color-scheme: dark)",
        "scheme": "slate",
        "primary": "deep-purple",
        "accent": "deep-purple",
        "toggle": {
            "icon": "material/brightness-4",
            "name": "Switch to light mode",
        }
    }],

    "site_url": "https://ipctk.xyz",

    # Set the repo location to get a badge with stats
    "repo_url": "https://github.com/ipc-sim/ipc-toolkit",
    "repo_name": "ipc-sim/ipc-toolkit",
    "icon": {"repo": "fontawesome/brands/github"},

    "edit_uri": "blob/main/docs/source",

    "features": [
        "navigation.expand",
        "navigation.tabs",
        "navigation.top",
        "navigation.tracking",
        "search.highlight",
        "search.share",
        "toc.follow",
        "content.tabs.link"
    ],

    "font": {
        "text": "Roboto",  # used for all the pages' text
        "code": "Roboto Mono"  # used for literal code blocks
    },

    "toc_title": "Contents",

    "version_dropdown": True,
    "version_json": "https://ipctk.xyz/versions.json",
}

html_title = "IPC Toolkit"

html_logo = "_static/hammer-wrench.svg"
html_favicon = "_static/favicon.ico"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = ["css/custom.css"]

# html_last_updated_fmt = "%B %d, %Y"
