Continuous Collision Detection
==============================

.. doxygenfunction:: ipc::is_step_collision_free

.. doxygenfunction:: ipc::compute_collision_free_stepsize

Narrow Phase CCD
----------------

.. doxygenclass:: ipc::NarrowPhaseCCD
    :allow-dot-graphs:

Tight Inclusion CCD
^^^^^^^^^^^^^^^^^^^

.. doxygenclass:: ipc::TightInclusionCCD
    :allow-dot-graphs:

Additive CCD
^^^^^^^^^^^^

.. doxygenclass:: ipc::AdditiveCCD
    :allow-dot-graphs:

Inexact CCD
^^^^^^^^^^^

.. note::
    This method is disabled by default. To enable it, set the
    ``IPC_TOOLKIT_WITH_INEXACT_CCD`` CMake option to ``ON``.

.. .. doxygenclass:: ipc::InexactCCD
..     :allow-dot-graphs:

.. doxygenfunction:: ipc::inexact_point_edge_ccd_2D

Nonlinear CCD
-------------

.. doxygenclass:: ipc::NonlinearTrajectory
    :allow-dot-graphs:

.. doxygenclass:: ipc::IntervalNonlinearTrajectory
    :allow-dot-graphs:


.. doxygenfunction:: ipc::point_point_nonlinear_ccd
.. doxygenfunction:: ipc::point_edge_nonlinear_ccd
.. doxygenfunction:: ipc::edge_edge_nonlinear_ccd
.. doxygenfunction:: ipc::point_triangle_nonlinear_ccd

Generic Interface
^^^^^^^^^^^^^^^^^

.. doxygenfunction:: ipc::conservative_piecewise_linear_ccd

Miscellaneous
-------------

.. doxygenfunction:: ipc::point_static_plane_ccd