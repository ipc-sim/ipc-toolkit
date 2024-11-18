Continuous Collision Detection
==============================

.. autofunction:: ipctk.is_step_collision_free

.. autofunction:: ipctk.compute_collision_free_stepsize

Narrow Phase CCD
----------------

.. autoclass:: ipctk.NarrowPhaseCCD

   .. autoclasstoc::

Tight Inclusion CCD
^^^^^^^^^^^^^^^^^^^

.. autoclass:: ipctk.TightInclusionCCD

   .. autoclasstoc::

Additive CCD
^^^^^^^^^^^^

.. autoclass:: ipctk.AdditiveCCD

   .. autoclasstoc::

Inexact CCD
^^^^^^^^^^^

.. note::
    This method is disabled by default. To enable it, set the
    ``IPC_TOOLKIT_WITH_INEXACT_CCD`` CMake option to ``ON``.

.. .. autoclass:: ipctk.InexactCCD

.. autofunction:: ipctk.inexact_point_edge_ccd_2D

Nonlinear CCD
-------------

.. autoclass:: ipctk.NonlinearTrajectory

   .. autoclasstoc::

.. autoclass:: ipctk.IntervalNonlinearTrajectory

    .. autoclasstoc::

.. autofunction:: ipctk.point_point_nonlinear_ccd
.. autofunction:: ipctk.point_edge_nonlinear_ccd
.. autofunction:: ipctk.edge_edge_nonlinear_ccd
.. autofunction:: ipctk.point_triangle_nonlinear_ccd

Miscellaneous
-------------

.. autofunction:: ipctk.point_static_plane_ccd