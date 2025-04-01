Barrier
=======

.. doxygenfunction:: barrier(const double, const double)
.. doxygenfunction:: barrier_first_derivative
.. doxygenfunction:: barrier_second_derivative

Barrier Force Magnitude
-----------------------

.. doxygenfunction:: barrier_force_magnitude
.. doxygenfunction:: barrier_force_magnitude_gradient

Adaptive Barrier Stiffness
--------------------------

.. doxygenfunction:: ipc::initial_barrier_stiffness
.. doxygenfunction:: ipc::update_barrier_stiffness

Semi-Implicit Stiffness
~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: ipc::semi_implicit_stiffness(const CollisionMesh&, const Eigen::MatrixXd&, const StencilsT&, const Eigen::VectorXd&, const Eigen::SparseMatrix<double>&, const double)

Barrier Class
-------------

.. doxygenclass:: ipc::Barrier
    :allow-dot-graphs:

Clamped Log Barrier
~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: ipc::ClampedLogBarrier
    :allow-dot-graphs:

Normalized Clamped Log Barrier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: ipc::NormalizedClampedLogBarrier
    :allow-dot-graphs:

Clamped Log Squared Barrier
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenclass:: ipc::ClampedLogSqBarrier
    :allow-dot-graphs:

Cubic Barrier
~~~~~~~~~~~~~

.. doxygenclass:: ipc::CubicBarrier
    :allow-dot-graphs: