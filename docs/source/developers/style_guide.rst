Style Guide
===========

This document provides a guide to the style used in this project.

Code Formatting
---------------

We utilize `ClangFormat <https://clang.llvm.org/docs/ClangFormat.html>`_ to automate code formatting. Please format your code before pushing and/or creating a pull request.

The project uses the root ``.clang-format`` (80 columns, WebKit-based).
Under ``tests/``, ``tests/.clang-format`` inherits that style and sets
``SortIncludes: false``. CI runs clang-format 20; format with the same version
locally to avoid formatting check failures (e.g. ``clang-format -i`` using
version 20, or use the pre-commit hook from :doc:`developers/tools`).
clang-tidy uses the same style via ``FormatStyle: file`` (see ``.clang-tidy``).

Additionally, ensure that your code adheres to the project's linting rules. Use the provided linting tools to check for any issues before committing your changes.

Naming conventions
------------------

General
^^^^^^^

In general, we stick to the following naming conventions:

* ``snake_case`` for variables, functions, and filenames
* ``PascalCase`` for classes and structs
* ``ALL_CAPS`` for constants and enum members
* ``m_`` prefix for class member variables (if the member is not public)
* ``member()`` to get a class member variable (if the member is not public)
* ``set_member()`` to set a class member variable (if the member is not public)

Specific
^^^^^^^^

* Vertex positions: ``vertices`` or ``positions``
* Vertex displacements: ``displacements``
* Vertex rest positions/material coordinates: ``rest_positions``
* Vertex velocities: ``velocities``
* Mesh edge matrix: ``edges``
* Mesh face matrix: ``faces``
* Element vertices: Use a numeral suffix (e.g., ``e0`` and ``e1`` for the end-points of an edge).
* Edge-edge pairings: Use suffixes ``a`` and ``b``.
* Continuous collision detection pairs: Use the suffix ``_t0`` for starting values and ``_t1`` for end values.
* Favor the term "collision" over "contact" (but this is not a hard rule).
* Prefer the term "potential" over "constraint" when referring to collisions and friction.

Documentation
-------------

We use `Doxygen <https://www.doxygen.nl/index.html>`_ to generate documentation. Please document your code before pushing and/or creating a pull request.