Barrier
=======

.. toctree::
   :caption: Barrier
   :hidden:

   Barrier.rst
   ClampedLogBarrier.rst
   NormalizedBarrier.rst
   ClampedLogSqBarrier.rst
   CubicBarrier.rst

Barrier functions and functionals.

Types
-----

.. list-table::
    :header-rows: 1

    * - Name
      - Description
    * - :doc:`Barrier <Barrier>`
      - Base class for barrier functions.
    * - :doc:`ClampedLogBarrier <ClampedLogBarrier>`
      - Smoothly clamped log barrier functions from [Li et al.].
    * - :doc:`ClampedLogSqBarrier <ClampedLogSqBarrier>`
      - Clamped log barrier with a quadratic log term from [Huang et al.].
    * - :doc:`CubicBarrier <CubicBarrier>`
      - Cubic barrier function from [Ando 2024].
    * - :doc:`NormalizedBarrier <NormalizedBarrier>`
      - Normalized barrier function from [Li et al.].

Functions
---------

.. list-table::
    :header-rows: 1

    * - Name
      - Description
    * - :func:`ipc::barrier`
      - Evaluate the barrier function.
    * - :func:`ipc::barrier_first_derivative`
      - Derivative of the barrier function.
    * - :func:`ipc::barrier_second_derivative`
      - Second derivative of the barrier function.
    * - :func:`ipc::barrier_force_magnitude`
      - Compute the barrier force magnitude.
    * - :func:`ipc::barrier_force_magnitude_gradient`
      - Compute the gradient of the barrier force magnitude.
    * - :func:`ipc::initial_barrier_stiffness`
      - Compute the initial barrier stiffness.
    * - :func:`ipc::update_barrier_stiffness`
      - Update the barrier stiffness based on the current state.
    * - :func:`ipc::semi_implicit_stiffness`
      - Compute the semi-implicit stiffness for all collisions.

Function Details
----------------

.. doxygenfunction:: barrier(const double, const double)
.. doxygenfunction:: barrier_first_derivative
.. doxygenfunction:: barrier_second_derivative

Barrier Force Magnitude
~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: barrier_force_magnitude
.. doxygenfunction:: barrier_force_magnitude_gradient

Adaptive Barrier Stiffness
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. doxygenfunction:: initial_barrier_stiffness
.. doxygenfunction:: update_barrier_stiffness

Semi-Implicit Stiffness
^^^^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: semi_implicit_stiffness(const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>, const StencilsT&, Eigen::ConstRef<Eigen::VectorXd>, const Eigen::SparseMatrix<double>&, const double)
