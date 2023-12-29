Style Guide
===========

This document provides a guide to the style used in this project.

Code Formatting
---------------

We utilize `ClangFormat <https://clang.llvm.org/docs/ClangFormat.html>`_ to automate code formatting. Please format your code before pushing and/or creating a pull request.

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

* vertex positions: ``vertices`` or ``positions``
* vertex displacements: ``displacements``
* vertex rest positions/material coordinates: ``rest_positions``
* vertex velocities: ``velocities``
* mesh edge matrix: ``edges``
* mesh face matrix: ``faces``
* element vertices: we use a numeral suffix (e.g., ``e0`` and ``e1`` for the end-points of an edge)
* edge-edge pairings: suffix of ``a`` and ``b``
* continuous collision detection pairs: suffix of ``_t0`` for starting values and ``_t1`` for end values
* favor the term "collision" over "contact" (but this is not a hard rule)
* we prefer the term "potential" over "constraint" when referring to the collisions and friction

Documentation
-------------

We use `Doxygen <https://www.doxygen.nl/index.html>`_ to generate documentation. Please document your code before pushing and/or creating a pull request.