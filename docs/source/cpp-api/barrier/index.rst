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
   TwoStageBarrier.rst

Barrier functions and functionals.

Types
-----

.. list-table::
    :header-rows: 1

    * - Name
      - Description
    * - :cpp:class:`ipc::Barrier`
      - Base class for barrier functions.
    * - :cpp:class:`ipc::ClampedLogBarrier`
      - Smoothly clamped log barrier functions from [Li et al.].
    * - :cpp:class:`ipc::ClampedLogSqBarrier`
      - Clamped log barrier with a quadratic log term from [Huang et al.].
    * - :cpp:class:`ipc::CubicBarrier`
      - Cubic barrier function from [Ando 2024].
    * - :cpp:class:`ipc::NormalizedBarrier`
      - Normalized barrier function from [Li et al.].
    * - :cpp:class:`ipc::TwoStageBarrier`
      - Two-stage barrier function from [Chen et al. 2025].


Functions
---------

.. list-table::
    :header-rows: 1

    * - Name
      - Description
    * - :cpp:func:`ipc::barrier`
      - Evaluate the barrier function.
    * - :cpp:func:`ipc::barrier_first_derivative`
      - Derivative of the barrier function.
    * - :cpp:func:`ipc::barrier_second_derivative`
      - Second derivative of the barrier function.
    * - :cpp:func:`ipc::barrier_force_magnitude`
      - Compute the barrier force magnitude.
    * - :cpp:func:`ipc::barrier_force_magnitude_gradient`
      - Compute the gradient of the barrier force magnitude.
    * - :cpp:func:`ipc::initial_barrier_stiffness`
      - Compute the initial barrier stiffness.
    * - :cpp:func:`ipc::update_barrier_stiffness`
      - Update the barrier stiffness based on the current state.
    * - :cpp:func:`ipc::semi_implicit_stiffness`
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
