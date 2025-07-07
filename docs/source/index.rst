.. toctree::
    :caption: About
    :hidden:

    Home <self>
    release_notes.rst
    license.rst
    dependencies.rst

.. toctree::
    :caption: Building
    :hidden:

    Library <building.rst>
    python.rst

.. toctree::
    :caption: Tutorial
    :hidden:

    tutorial/getting_started.rst
    tutorial/convergent.rst
    tutorial/nonlinear_ccd.rst
    tutorial/adhesion.rst
    tutorial/simulation.rst
    tutorial/misc.rst
    tutorial/faq.rst
    tutorial/references.rst

.. toctree::
    :caption: C++ API
    :hidden:

    cpp-api/potentials.rst
    cpp-api/collision_mesh.rst
    cpp-api/candidates.rst
    cpp-api/normal_collisions.rst
    cpp-api/tangential_collisions.rst
    cpp-api/friction.rst
    cpp-api/broad_phase.rst
    cpp-api/ccd.rst
    cpp-api/distance.rst
    cpp-api/tangent.rst
    cpp-api/barrier.rst
    cpp-api/adhesion.rst
    cpp-api/intersections.rst
    cpp-api/interval.rst
    cpp-api/utils.rst

.. toctree::
    :caption: Python API
    :hidden:

    python-api/potentials.rst
    python-api/collision_mesh.rst
    python-api/candidates.rst
    python-api/normal_collisions.rst
    python-api/tangential_collisions.rst
    python-api/friction.rst
    python-api/broad_phase.rst
    python-api/ccd.rst
    python-api/distance.rst
    python-api/tangent.rst
    python-api/barrier.rst
    python-api/adhesion.rst
    python-api/intersections.rst
    python-api/interval.rst
    python-api/utils.rst

.. toctree::
    :caption: Developers
    :hidden:

    contributing
    style_guide
    Code of Conduct <code_of_conduct.md>

.. image:: _static/logo.png
    :alt: IPC Toolkit

.. raw:: html

    <div align="center">

.. image:: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml/badge.svg
    :target: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml
    :alt: Continuous Integration
.. image:: https://codecov.io/github/ipc-sim/ipc-toolkit/graph/badge.svg?token=9BR6GPKRY8
    :target: https://codecov.io/github/ipc-sim/ipc-toolkit
    :alt: Coverage
.. image:: https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue
    :target: https://github.com/ipc-sim/ipc-toolkit/blob/main/LICENSE
    :alt: License
.. image:: https://img.shields.io/pypi/v/ipctk?color=brightgreen&label=PyPI&logo=python&logoColor=white
    :target: https://pypi.org/project/ipctk/
    :alt: PyPI
.. image:: https://img.shields.io/pypi/dm/ipctk?label=PyPI%20Downloads&logo=python&logoColor=white
    :target: https://pypi.org/project/ipctk/
    :alt: PyPI Downloads

.. raw:: html

    </div>

IPC Toolkit is a set of reusable functions to integrate Incremental Potential Contact (IPC) into a simulation.

**Features**

- IPC barrier function and its derivatives and adaptive barrier stiffness algorithm
- Broad- and narrow-phase continuous collision detection (CCD) of linear and nonlinear trajectories
- Distance computation and derivatives between edges in 2D and triangles in 3D
- Distance barrier potential and its derivatives
- Smooth and lagged dissipative friction potential and its derivatives

**Limitations**

This is not a full simulation library. As such it does not include any physics or solvers. For a full simulation implementation, we recommend `PolyFEM <https://polyfem.github.io/>`_ (a finite element library) or `Rigid IPC <https://github.com/ipc-sim/rigid-ipc>`_ (rigid-body dynamics) both of which utilize the IPC Toolkit.

**Usage**

See the `tutorial <https://ipctk.xyz/tutorial/getting_started.html>`_ for a quick introduction to the toolkit, or the `documentation <https://ipctk.xyz/cpp.html>`_ for a full reference.

**Python Bindings**

We provide Python bindings for functions in the toolkit using `pybind11 <https://github.com/pybind/pybind11>`_. For more information see the `Python documentation <https://ipctk.xyz/python.html>`_.

**Contributing**

This project is open to contributors! Contributions can come in the form of feature requests, bug fixes, documentation, tutorials, and the like. We highly recommend filing an Issue first before submitting a Pull Request.

Simply fork this repository and make a Pull Request! We would appreciate:

* Implementation of new features
* Bug Reports
* Documentation
* Testing

**Citation**

IPC Toolkit is created and maintained by academics: citations let us know our work is having impact! Please cite the IPC Toolkit or otherwise give a shout-out if and when it contributes to published works.

.. code-block:: bibtex

    @software{ipc_toolkit,
        author = {Zachary Ferguson and others},
        title = {{IPC Toolkit}},
        url = {https://github.com/ipc-sim/ipc-toolkit},
        year = {2020},
    }

**License**

MIT License © 2020, the IPC-Sim organization (See `License <license.html>`_ for details).