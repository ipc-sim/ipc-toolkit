Continuous Collision Detection
==============================

.. doxygenfunction:: ipc::is_step_collision_free

.. doxygenfunction:: ipc::compute_collision_free_stepsize

.. doxygenvariable:: ipc::DEFAULT_CCD_TOLERANCE
.. doxygenvariable:: ipc::DEFAULT_CCD_MAX_ITERATIONS
.. doxygenvariable:: ipc::DEFAULT_CCD_CONSERVATIVE_RESCALING

Individual CCD Functions
------------------------

.. doxygenfunction:: ipc::point_point_ccd
.. doxygenfunction:: ipc::point_edge_ccd
.. doxygenfunction:: ipc::edge_edge_ccd
.. doxygenfunction:: ipc::point_triangle_ccd

Generic Interface
^^^^^^^^^^^^^^^^^

.. doxygenfunction:: ipc::ccd_strategy

Additive CCD
------------

.. doxygenfunction:: ipc::additive_ccd::point_point_ccd
.. doxygenfunction:: ipc::additive_ccd::point_edge_ccd
.. doxygenfunction:: ipc::additive_ccd::edge_edge_ccd
.. doxygenfunction:: ipc::additive_ccd::point_triangle_ccd

Generic Interface
^^^^^^^^^^^^^^^^^

.. doxygenfunction:: ipc::additive_ccd::additive_ccd

Nonlinear CCD
-------------

.. doxygenclass:: ipc::NonlinearTrajectory
.. doxygenclass:: ipc::IntervalNonlinearTrajectory

.. doxygenfunction:: ipc::point_point_nonlinear_ccd
.. doxygenfunction:: ipc::point_edge_nonlinear_ccd
.. doxygenfunction:: ipc::edge_edge_nonlinear_ccd
.. doxygenfunction:: ipc::point_triangle_nonlinear_ccd

Generic Interface
^^^^^^^^^^^^^^^^^

.. doxygenfunction:: ipc::conservative_piecewise_linear_ccd