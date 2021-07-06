Incremental Potential Contact (IPC) Toolkit
===========================================
.. image:: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml/badge.svg
   :target: https://github.com/ipc-sim/ipc-toolkit/actions/workflows/continuous.yml
   :alt: Build
.. image:: https://img.shields.io/github/license/ipc-sim/ipc-toolkit.svg?color=blue
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT

A set of reusable functions to integrate IPC into an existing simulation.

- **Free Software**: MIT License
- **Github Repository**: https://github.com/ipc-sim/ipc-toolkit
- **Paper**: https://ipc-sim.github.io/file/IPC-paper-fullRes.pdf

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: General

    .. installation
    .. concepts
    .. examples
    changelog

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: Developers

    contributing

.. toctree::
    :hidden:
    :maxdepth: 2
    :caption: API Documentation

    api/library_root

API
---

Barrier Function
~~~~~~~~~~~~~~~~

:func:`ipc::barrier`

:func:`ipc::barrier_gradient`

:func:`ipc::barrier_hessian`

:func:`ipc::point_edge_distance_type`

Development Status
------------------

- ✅ contacts
- ✅ friction

.. warning::
    This toolkit is in an early stage of development and consequently the API may
    change to better improve usability. If you have any problems
    or find any bugs please `post and issue <https://github.com/ipc-sim/ipc-toolkit/issues>`__. For a complete list of changes,
    please see `changelog <changelog.html>`__.
    Meanwhile, for a definitive reference for these functions, please see the `IPC source code <https://github.com/ipc-sim/IPC>`__.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
