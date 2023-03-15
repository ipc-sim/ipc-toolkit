Friction
========

Friction Constraints
--------------------

.. doxygenclass:: ipc::FrictionConstraints

Friction Constraint
^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: ipc::FrictionConstraint

Vertex-Vertex Friction Constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: ipc::VertexVertexFrictionConstraint

Edge-Vertex Friction Constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: ipc::EdgeVertexFrictionConstraint

Edge-Edge Friction Constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: ipc::EdgeEdgeFrictionConstraint

Triangle-Vertex Friction Constraint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: ipc::FaceVertexFrictionConstraint

Smooth Mollifier
----------------

.. doxygenfunction:: f0_SF
.. doxygenfunction:: f1_SF_over_x
.. doxygenfunction:: df1_x_minus_f1_over_x3

Normal Force Magnitude
----------------------

.. doxygenfunction:: ipc::compute_normal_force_magnitude
.. doxygenfunction:: ipc::compute_normal_force_magnitude_gradient

Tangent Basis
-------------

.. doxygenfunction:: ipc::point_point_tangent_basis
.. doxygenfunction:: ipc::point_edge_tangent_basis
.. doxygenfunction:: ipc::edge_edge_tangent_basis
.. doxygenfunction:: ipc::point_triangle_tangent_basis

Tangent Basis Jacobians
^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: ipc::point_point_tangent_basis_jacobian
.. doxygenfunction:: ipc::point_edge_tangent_basis_jacobian
.. doxygenfunction:: ipc::edge_edge_tangent_basis_jacobian
.. doxygenfunction:: ipc::point_triangle_tangent_basis_jacobian

Relative Velocity
-----------------

.. doxygenfunction:: ipc::point_point_relative_velocity
.. doxygenfunction:: ipc::point_edge_relative_velocity
.. doxygenfunction:: ipc::edge_edge_relative_velocity
.. doxygenfunction:: ipc::point_triangle_relative_velocity

Relative Velocity as Multiplier Matricies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: ipc::point_point_relative_velocity_matrix
.. doxygenfunction:: ipc::point_edge_relative_velocity_matrix
.. doxygenfunction:: ipc::edge_edge_relative_velocity_matrix
.. doxygenfunction:: ipc::point_triangle_relative_velocity_matrix

Relative Velocity Matrix Jacobians
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: ipc::point_point_relative_velocity_matrix_jacobian
.. doxygenfunction:: ipc::point_edge_relative_velocity_matrix_jacobian
.. doxygenfunction:: ipc::edge_edge_relative_velocity_matrix_jacobian
.. doxygenfunction:: ipc::point_triangle_relative_velocity_matrix_jacobian