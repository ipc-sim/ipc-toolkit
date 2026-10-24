# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import contextlib
import io
import logging
import os
import pathlib
import subprocess
import sys

# A full build narrates itself for a few thousand lines, which buries the
# handful of warnings worth acting on. We keep it quiet by default and let
# DOCS_VERBOSE=1 bring the whole log back when a build needs debugging.
VERBOSE = bool(os.environ.get("DOCS_VERBOSE"))


def quiet_stdout():
    """Swallow stdout from code that writes to it directly.

    Sphinx's own -q only silences messages that go through its status stream,
    so anything using a bare print() needs to be wrapped by hand.
    """
    if VERBOSE:
        return contextlib.nullcontext()
    return contextlib.redirect_stdout(io.StringIO())


# -- Progress reporting ------------------------------------------------------
# A quiet build says nothing for the half minute or so it runs, which looks
# exactly like a hung one. So each slow stage spins a single self-erasing line
# naming what it is doing, and stamps that line with a result when it finishes.
#
# This only happens on an interactive terminal. Redirected output (CI, a log
# file) stays as quiet as it was, and a verbose build skips the spinner because
# its own logs already show progress.
try:
    from yaspin import yaspin
except ImportError:  # the docs still build without it, just without progress
    yaspin = None

SHOW_PROGRESS = yaspin is not None and not VERBOSE and sys.stdout.isatty()


class SpinnerSafeStream:
    """Keep the spinner from scribbling over whatever is written to a stream.

    The spinner owns the last line of the terminal and redraws it on a timer,
    so a warning written underneath lands in the middle of that line. We clear
    the spinner first, write, then let it resume on a fresh line.
    """

    def __init__(self, stream, spinner):
        self._stream = stream
        self._spinner = spinner

    def write(self, text):
        try:
            self._spinner.hide()
            written = self._stream.write(text)
            self._stream.flush()
        finally:
            self._spinner.show()
        return written

    def __getattr__(self, name):
        return getattr(self._stream, name)


@contextlib.contextmanager
def stage(text):
    """Spin while one build stage runs, then stamp the line with its result."""
    if not SHOW_PROGRESS:
        yield None
        return
    with yaspin(text=text, color="yellow") as spinner:
        try:
            yield spinner
        except BaseException:
            spinner.fail("\U0001f4a5 ")
            raise
        spinner.ok("\u2705 ")


# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
# -- Project information -----------------------------------------------------
from datetime import datetime

sys.path.append(str(pathlib.Path(__file__).parents[2] / "python"))
with quiet_stdout():
    from _find_ipctk import ipctk  # noqa

project = "IPC Toolkit"
copyright = f"2020-{datetime.now().year}, Zachary Ferguson; MIT License"
author = "Zachary Ferguson"
version = ipctk.__version__

# -- General configuration ---------------------------------------------------

# Doxygen
# -q silences Doxygen's per-file progress. Its warnings go to stderr and are
# unaffected, so we still hear about anything that is actually wrong. While the
# spinner is up we hold those warnings in a pipe and replay them once the
# spinner has cleared its line, so they arrive intact rather than interleaved.
pathlib.Path("../build/doxyoutput").mkdir(parents=True, exist_ok=True)
doxygen_flags = [] if VERBOSE else ["-q"]

with stage("Doxygen") as doxygen_spinner:
    doxygen = subprocess.run(
        ["doxygen", *doxygen_flags, "Doxyfile"],
        stderr=subprocess.PIPE if doxygen_spinner is not None else None,
    )
    if doxygen.stderr:
        doxygen_spinner.hide()
        sys.stderr.buffer.write(doxygen.stderr)
        sys.stderr.flush()
        doxygen_spinner.show()
    if doxygen.returncode != 0:
        raise SystemExit("Doxygen failed! Exiting")

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

# sphinx_immaterial announces where it wrote the sitemap with a bare print(),
# so -q never reaches it. We swap in a wrapper around the real handler here,
# at config-read time, because that happens before Sphinx loads the extension
# and binds the handler by name. If the theme ever renames the function the
# getattr below just returns None, and the worst case is that one line comes
# back rather than the build failing.
if not VERBOSE:
    try:
        from sphinx_immaterial import postprocess_html
    except ImportError:
        postprocess_html = None

    create_sitemap = getattr(postprocess_html, "create_sitemap", None)
    if create_sitemap is not None:

        def quiet_create_sitemap(app, exception, _wrapped=create_sitemap):
            with quiet_stdout():
                _wrapped(app, exception)

        postprocess_html.create_sitemap = quiet_create_sitemap

bibtex_bibfiles = ["references.bib"]
bibtex_reference_style = "author_year"
bibtex_default_style = "plain"

myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
]

object_description_options = [
    ("cpp:.*", dict(clang_format_style={"BasedOnStyle": "WebKit"})),
]

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

root_doc = "index"

suppress_warnings = ["myst.header", "autodoc.import_object", "autosummary"]

# Setup the breathe extension
breathe_projects = {project: "../build/doxyoutput/xml"}
breathe_default_project = project
breathe_default_members = (
    "members",
    "undoc-members",
    "protected-members",
    "private-members",
)
breathe_show_define_initializer = True
breathe_show_include = True

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": True,
    "special-members": True,
    "show-inheritance": True,
}

# -- GraphViz configuration ----------------------------------
graphviz_output_format = "svg"

graphviz_dot_args = [
    "-Ecolor=#CE93D8",
    "-Kdot",
    "-Gbgcolor=transparent",
    "-Nfontname=Menlo",
]

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
    "palette": [
        {
            "media": "(prefers-color-scheme: light)",
            "scheme": "default",
            "primary": "deep-purple",
            "accent": "deep-purple",
            "toggle": {
                "icon": "material/brightness-7",
                "name": "Switch to dark mode",
            },
        },
        {
            "media": "(prefers-color-scheme: dark)",
            "scheme": "slate",
            "primary": "deep-purple",
            "accent": "deep-purple",
            "toggle": {
                "icon": "material/brightness-4",
                "name": "Switch to light mode",
            },
        },
    ],
    "site_url": "https://ipctk.xyz",
    # Set the repo location to get a badge with stats
    "repo_url": "https://github.com/ipc-sim/ipc-toolkit",
    "repo_name": "ipc-sim/ipc-toolkit",
    "icon": {"repo": "fontawesome/brands/github"},
    "features": [
        "content.tabs.link",
        "navigation.footer",
        "navigation.tabs",
        "navigation.top",
        "navigation.tracking",
        "search.highlight",
        "search.share",
        "toc.follow",
    ],
    "font": {
        "text": "Roboto",  # used for all the pages' text
        "code": "Roboto Mono",  # used for literal code blocks
    },
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

# -- Progress reporting for Sphinx's own stages -------------------------------


def connect_progress(app):
    """Drive one spinner through Sphinx's read and write passes.

    Those two passes are the slow half of the build and -q makes them silent,
    so we name the running phase and count the documents as they go by.

    Warnings have to keep printing cleanly underneath all that. Sphinx binds
    its warning stream when the application is constructed, which is before any
    event we can hook, so swapping sys.stderr is not enough on its own -- we
    also wrap the stream on Sphinx's own warning handler. Both wrappers clear
    the spinner before writing, so a warning never lands mid-frame.
    """
    state = {
        "spinner": None,
        "stderr": None,
        "handlers": [],
        "phase": "",
        "done": 0,
        "total": 0,
    }

    def retitle():
        spinner = state["spinner"]
        if spinner is None:
            return
        progress = ""
        if state["total"]:
            progress = f" ({state['done']}/{state['total']})"
        elif state["done"]:
            progress = f" ({state['done']})"
        spinner.text = f"Sphinx{state['phase']}{progress}"

    def warning_handlers():
        """Sphinx's warning handler, the one thing above the WARNING level."""
        return [
            handler
            for handler in logging.getLogger("sphinx").handlers
            if handler.level >= logging.WARNING and hasattr(handler, "stream")
        ]

    def on_start(_app):
        spinner = yaspin(text="Sphinx", color="yellow")
        spinner.start()
        state["spinner"] = spinner
        state["stderr"] = sys.stderr
        sys.stderr = SpinnerSafeStream(sys.stderr, spinner)
        for handler in warning_handlers():
            state["handlers"].append((handler, handler.stream))
            handler.stream = SpinnerSafeStream(handler.stream, spinner)

    def on_read_start(_app, _env, docnames):
        state.update(phase=": reading sources", done=0, total=len(docnames))
        retitle()

    def on_read_doc(_app, _docname, _source):
        state["done"] += 1
        retitle()

    def on_write_page(_app, _pagename, _templatename, _context, _doctree):
        if not state["phase"].endswith("writing output"):
            state.update(phase=": writing output", done=0, total=0)
        state["done"] += 1
        retitle()

    def on_finish(_app, exception):
        spinner = state["spinner"]
        if spinner is None:
            return
        state["spinner"] = None
        sys.stderr = state["stderr"]
        for handler, stream in state["handlers"]:
            handler.stream = stream
        state["handlers"].clear()
        spinner.text = "Sphinx"
        if exception is None:
            spinner.ok("\u2705 ")
        else:
            spinner.fail("\U0001f4a5 ")

    app.connect("builder-inited", on_start)
    app.connect("env-before-read-docs", on_read_start)
    app.connect("source-read", on_read_doc)
    app.connect("html-page-context", on_write_page)
    # Last, so the spinner covers everything else that runs at the end.
    app.connect("build-finished", on_finish, priority=900)


# -- Custom skip logic for autodoc --------------------------------------------


def setup(app):
    if SHOW_PROGRESS:
        connect_progress(app)

    def skip(app, what, name, obj, skip, options):
        # Skip the specific private attribute that is causing the crash
        if name == "__entries":
            return True
        # Skip the internal pybind11 base classes and builtins
        if name == "pybind11_object" or name.startswith("pybind11_builtins"):
            return True
        # Optional: Skip the module's "self" capsule if it appears
        if name == "PyCapsule":
            return True
        return skip

    # Connect the function to the event
    app.connect("autodoc-skip-member", skip)
