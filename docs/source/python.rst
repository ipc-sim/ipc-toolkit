Build and Install
=================

We provide Python bindings for functions in the toolkit using
`pybind11`_.

.. _pybind11: https://github.com/pybind/pybind11

Currently, the bindings must be built from scratch. The easiest way to
do this is to use the ``setup.py`` script which uses ``setuptools``. For
example:

.. code:: sh

   python setup.py install

will build the library and python bindings and then install them on your
system.

You can test the install was successful by doing
``python -c "import ipctk"``.

API
===

Main Functions
--------------

.. autofunction:: ipctk.construct_constraint_set

.. autofunction:: ipctk.compute_barrier_potential
.. autofunction:: ipctk.compute_barrier_potential_gradient
.. autofunction:: ipctk.compute_barrier_potential_hessian

Collision Constraints
---------------------

.. autoclass:: ipctk.CollisionConstraint
    :members:
    :show-inheritance:
.. autoclass:: ipctk.VertexVertexConstraint
    :members:
    :show-inheritance:
.. autoclass:: ipctk.EdgeVertexConstraint
    :members:
    :show-inheritance:
.. autoclass:: ipctk.EdgeEdgeConstraint
    :members:
    :show-inheritance:
.. autoclass:: ipctk.FaceVertexConstraint
    :members:
    :show-inheritance:
.. autoclass:: ipctk.Constraints
    :members:

Barrier
-------

.. autofunction:: ipctk.barrier

.. autofunction:: ipctk.barrier_gradient

.. autofunction:: ipctk.barrier_hessian

Distance
--------

Distance Type
^^^^^^^^^^^^^

.. autoclass:: ipctk.PointEdgeDistanceType
.. autoclass:: ipctk.EdgeEdgeDistanceType
.. autoclass:: ipctk.PointTriangleDistanceType

.. autofunction:: ipctk.point_edge_distance_type
.. autofunction:: ipctk.edge_edge_distance_type
.. autofunction:: ipctk.point_triangle_distance_type

Edge-Edge Mollifier
^^^^^^^^^^^^^^^^^^^

.. autofunction:: ipctk.edge_edge_mollifier_threshold
.. autofunction:: ipctk.edge_edge_cross_squarednorm
.. autofunction:: ipctk.edge_edge_cross_squarednorm_gradient
.. autofunction:: ipctk.edge_edge_cross_squarednorm_hessian
.. autofunction:: ipctk.edge_edge_mollifier
.. autofunction:: ipctk.edge_edge_mollifier_gradient
.. autofunction:: ipctk.edge_edge_mollifier_hessian

Edge-Edge
^^^^^^^^^

.. autofunction:: ipctk.edge_edge_distance
.. autofunction:: ipctk.edge_edge_distance_gradient
.. autofunction:: ipctk.edge_edge_distance_hessian

Line-Line
^^^^^^^^^

.. autofunction:: ipctk.line_line_distance
.. autofunction:: ipctk.line_line_distance_gradient
.. autofunction:: ipctk.line_line_distance_hessian

Point-Edge
^^^^^^^^^^^

.. autofunction:: ipctk.point_edge_distance
.. autofunction:: ipctk.point_edge_distance_gradient
.. autofunction:: ipctk.point_edge_distance_hessian

Point-Line
^^^^^^^^^^^

.. autofunction:: ipctk.point_line_distance
.. autofunction:: ipctk.point_line_distance_gradient
.. autofunction:: ipctk.point_line_distance_hessian

Point-Plane
^^^^^^^^^^^

.. autofunction:: ipctk.point_plane_distance
.. autofunction:: ipctk.point_plane_distance_gradient
.. autofunction:: ipctk.point_plane_distance_hessian

Point-Point
^^^^^^^^^^^

.. autofunction:: ipctk.point_point_distance
.. autofunction:: ipctk.point_point_distance_gradient
.. autofunction:: ipctk.point_point_distance_hessian

Point-Triangle
^^^^^^^^^^^^^^

.. autofunction:: ipctk.point_triangle_distance
.. autofunction:: ipctk.point_triangle_distance_gradient
.. autofunction:: ipctk.point_triangle_distance_hessian

CCD
---

Broad-Phase
^^^^^^^^^^^

.. autoclass:: ipctk.BroadPhaseMethod

Candidates
""""""""""

.. autoclass:: ipctk.VertexVertexCandidate
    :members:
.. autoclass:: ipctk.EdgeVertexCandidate
    :members:
.. autoclass:: ipctk.EdgeEdgeCandidate
    :members:
.. autoclass:: ipctk.EdgeFaceCandidate
    :members:
.. autoclass:: ipctk.FaceVertexCandidate
    :members:

.. autoclass:: ipctk.Candidates

Narrow-Phase
^^^^^^^^^^^^

Utils
-----

.. autofunction:: ipctk.has_intersections

Logger
^^^^^^

.. autoclass:: ipctk.LoggerLevel
.. autofunction:: ipctk.set_logger_level

Multi-Threading
^^^^^^^^^^^^^^^

.. autofunction:: ipctk.get_num_threads
.. autofunction:: ipctk.set_num_threads

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`
