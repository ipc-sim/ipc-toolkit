# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# A full build prints roughly three thousand lines of Doxygen and Sphinx
# progress, which buries the handful of warnings we actually care about. So we
# build quietly by default: -q drops Sphinx's progress log while leaving
# warnings and errors on stderr, and conf.py passes the same flag to Doxygen.
# Set DOCS_VERBOSE=1 when a build misbehaves and you want the play-by-play:
#
#     make html DOCS_VERBOSE=1
#
ifeq ($(DOCS_VERBOSE),)
QUIET         = -q
endif

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile clean

clean:
	rm -rf source/api/
	@$(SPHINXBUILD) -M clean "$(SOURCEDIR)" "$(BUILDDIR)" $(QUIET) $(SPHINXOPTS) $(O)

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(QUIET) $(SPHINXOPTS) $(O)

xml:
	cd "$(SOURCEDIR)" && doxygen $(QUIET) Doxyfile

livehtml:
	sphinx-autobuild --host 0.0.0.0 -j auto -a -q "$(SOURCEDIR)" "$(BUILDDIR)/html" $(SPHINXOPTS) $(O)
