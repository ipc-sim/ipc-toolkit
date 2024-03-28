Continuous Collision Detection
==============================

.. autofunction:: ipctk.is_step_collision_free

.. autofunction:: ipctk.compute_collision_free_stepsize

Individual CCD Functions
------------------------

.. .. autovariable:: ipctk.DEFAULT_CCD_TOLERANCE
.. .. autovariable:: ipctk.DEFAULT_CCD_MAX_ITERATIONS
.. .. autovariable:: ipctk.DEFAULT_CCD_CONSERVATIVE_RESCALING

.. autofunction:: ipctk.point_point_ccd
.. autofunction:: ipctk.point_edge_ccd
.. autofunction:: ipctk.edge_edge_ccd
.. autofunction:: ipctk.point_triangle_ccd

Generic Interface
^^^^^^^^^^^^^^^^^

.. autofunction:: ipctk.ccd_strategy

Tight Inclusion CCD
-------------------

.. autofunction:: ipctk.tight_inclusion.edge_edge_ccd
.. autofunction:: ipctk.tight_inclusion.point_triangle_ccd

.. autofunction:: ipctk.tight_inclusion.compute_ccd_filters

Additive CCD
------------

.. autofunction:: ipctk.additive_ccd.point_point_ccd
.. autofunction:: ipctk.additive_ccd.point_edge_ccd
.. autofunction:: ipctk.additive_ccd.edge_edge_ccd
.. autofunction:: ipctk.additive_ccd.point_triangle_ccd

Generic Interface
^^^^^^^^^^^^^^^^^

.. autofunction:: ipctk.additive_ccd.additive_ccd

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

Generic Interface
^^^^^^^^^^^^^^^^^

.. autofunction:: ipctk.conservative_piecewise_linear_ccd