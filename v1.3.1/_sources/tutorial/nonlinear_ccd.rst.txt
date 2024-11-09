Nonlinear CCD
=============

We can also perform CCD of nonlinear trajectories (of linear geometry) using the method of :cite:t:`Ferguson2021RigidIPC`. While :cite:t:`Ferguson2021RigidIPC` introduce their method in the context of rigid bodies, it can be applied to any nonlinear trajectory of linear geometry. The method works by transforming the nonlinear trajectories into (adaptive) piecewise linear trajectories with an envelope/minimum separation around each piece, enclosing the nonlinear trajectory. The method then performs CCD on the piecewise linear trajectories to find the earliest time of impact.

We provide the following functions to perform nonlinear CCD:

.. md-tab-set::

    .. md-tab-item:: C++

        * :cpp:func:`ipc::point_point_nonlinear_ccd`,
        * :cpp:func:`ipc::point_edge_nonlinear_ccd`,
        * :cpp:func:`ipc::edge_edge_nonlinear_ccd`, and
        * :cpp:func:`ipc::point_triangle_nonlinear_ccd`.

        Each of these functions take as input a :cpp:func:`ipc::NonlinearTrajectory` object for the endpoints of the linear geometry.

    .. md-tab-item:: Python

        * :py:func:`ipctk.point_point_nonlinear_ccd`,
        * :py:func:`ipctk.point_edge_nonlinear_ccd`,
        * :py:func:`ipctk.edge_edge_nonlinear_ccd`, and
        * :py:func:`ipctk.point_triangle_nonlinear_ccd`.

        Each of these functions take as input a :py:func:`ipctk.NonlinearTrajectory` object for the endpoints of the linear geometry.

For example, the following code defines a rigid trajectory in 2D in order to perform nonlinear CCD between a point and edge:

.. md-tab-set::

    .. md-tab-item:: C++

        .. literalinclude:: ../../../tests/src/tests/ccd/test_nonlinear_ccd.cpp
            :language: c++
            :start-after: // BEGIN_RIGID_2D_TRAJECTORY
            :end-before: // END_RIGID_2D_TRAJECTORY

    .. md-tab-item:: Python

        .. literalinclude:: ../../../python/tests/test_ccd.py
            :language: python
            :start-after: # BEGIN_RIGID_2D_TRAJECTORY
            :end-before: # END_RIGID_2D_TRAJECTORY
            :dedent: 4

Defining the Trajectory
-----------------------

Let's dive deeper by breaking down the implementation of Rigid2DTrajectory. The first function we need to implement is the call operator:

.. md-tab-set::

    .. md-tab-item:: C++

        .. literalinclude:: ../../../tests/src/tests/ccd/test_nonlinear_ccd.cpp
            :language: c++
            :lines: 127-134
            :dedent: 4


    .. md-tab-item:: Python

        .. literalinclude:: ../../../python/tests/test_ccd.py
            :language: python
            :lines: 91-95
            :dedent: 8

This function computes the position of the point at a time :math:`t \in [0, 1]`. This defines the trajectory of the point. In this case, we have a rigid body with a center of mass (COM) at the origin. The trajectory of the point is given by:

.. math::

    x(t) := R(\theta + t \Delta \theta) \bar{x} + T + t \Delta T

where :math:`\theta` is the angle of rotation about the COM, :math:`T` is the translation of the COM, :math:`\Delta \theta` and :math:`\Delta T` are the updates to :math:`\theta` and :math:`T`, respectively, and :math:`\bar{x}` is the position of the point in the interial frame.

Computing a Conservative Envelope
---------------------------------

The second function we need to implement is ``max_distance_from_linear``.

.. md-tab-set::

    .. md-tab-item:: C++

        .. literalinclude:: ../../../tests/src/tests/ccd/test_nonlinear_ccd.cpp
            :language: c++
            :lines: 136-147
            :dedent: 4

    .. md-tab-item:: Python

        .. literalinclude:: ../../../python/tests/test_ccd.py
            :language: python
            :lines: 97-103
            :dedent: 8

This function computes the maximum distance over a time interval :math:`[t_0, t_1]` between the nonlinear trajectory and a line segment from :math:`x(t_0)` to :math:`x(t_1)`. Mathematically this function computes

.. math::

    \min_{t\in[0, 1]} \|x((t_1 - t_0) t + t_0) - ((x(t_1) - x(t_0))t + x(t_0))\|,

for a given start and end time :math:`t_0` and :math:`t_1`, respectively.

In the case of a 2D rigid body, we can compute this value analytically because we know the :math:`\arg\!\min`:

.. math::

    \underset{t\in[0, 1]}{\arg\!\min} \|x((t_1 - t_0) t + t_0) - ((x(t_1) - x(t_0))t + x(t_0))\| = 0.5,

for :math:`(t_1 - t_0) \Delta \theta \leq \pi/2`, otherwise we can use the most conservative envelope radius of :math:`2 \|\bar{x}\|`.

Performing Nonlinear CCD
------------------------

Last, we use the ``Rigid2DTrajectory`` to perform nonlinear CCD between a point and edge:

.. md-tab-set::

    .. md-tab-item:: C++

        .. literalinclude:: ../../../tests/src/tests/ccd/test_nonlinear_ccd.cpp
            :language: c++
            :start-after: // BEGIN_TEST_RIGID_2D_TRAJECTORY
            :end-before: // END_TEST_RIGID_2D_TRAJECTORY
            :dedent: 4

    .. md-tab-item:: Python

        .. literalinclude:: ../../../python/tests/test_ccd.py
            :language: python
            :start-after: # BEGIN_TEST_RIGID_2D_TRAJECTORY
            :end-before: # END_TEST_RIGID_2D_TRAJECTORY
            :dedent: 4

.. note::
    We adjust the ``conservative_rescaling`` parameter to get a more accurate time of impact (TOI), but in practice, this is not needed as a more conservative estimate of the TOI is sufficient to avoid penetrations.