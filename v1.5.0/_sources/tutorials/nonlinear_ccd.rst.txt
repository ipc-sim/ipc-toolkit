Nonlinear CCD
=============

We also implement CCD of nonlinear trajectories (of linear geometry) using the method of :cite:t:`Ferguson2021RigidIPC`. While :cite:t:`Ferguson2021RigidIPC` introduce their method in the context of rigid bodies, it can be applied to any nonlinear trajectory of linear geometry.

The method works by transforming the nonlinear trajectories into (adaptive) piecewise linear trajectories with an envelope/minimum separation around each piece, enclosing the nonlinear trajectory. The method then performs CCD on the piecewise linear trajectories to find the earliest time of impact.

We provide the following functions to perform nonlinear CCD:

.. md-tab-set::

    .. md-tab-item:: C++

        * :cpp:func:`ipc::point_point_nonlinear_ccd`,
        * :cpp:func:`ipc::point_edge_nonlinear_ccd`,
        * :cpp:func:`ipc::edge_edge_nonlinear_ccd`, and
        * :cpp:func:`ipc::point_triangle_nonlinear_ccd`.

        Each of these functions take as input a :cpp:class:`ipc::NonlinearTrajectory` object for the endpoints of the linear geometry.

    .. md-tab-item:: Python

        * :py:func:`ipctk.point_point_nonlinear_ccd`,
        * :py:func:`ipctk.point_edge_nonlinear_ccd`,
        * :py:func:`ipctk.edge_edge_nonlinear_ccd`, and
        * :py:func:`ipctk.point_triangle_nonlinear_ccd`.

        Each of these functions take as input a :py:class:`ipctk.NonlinearTrajectory` object for the endpoints of the linear geometry.

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

Let's dive deeper by breaking down the implementation of ``Rigid2DTrajectory``. The first function we need to implement is the call operator:

.. md-tab-set::

    .. md-tab-item:: C++

        .. literalinclude:: ../../../tests/src/tests/ccd/test_nonlinear_ccd.cpp
            :language: c++
            :start-after: // BEGIN_RIGID_2D_CALL
            :end-before: // END_RIGID_2D_CALL
            :dedent: 4


    .. md-tab-item:: Python

        .. literalinclude:: ../../../python/tests/test_ccd.py
            :language: python
            :start-after: # BEGIN_RIGID_2D_CALL
            :end-before: # END_RIGID_2D_CALL
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
            :start-after: // BEGIN_RIGID_2D_MAX_DISTANCE_FROM_LINEAR
            :end-before: // END_RIGID_2D_MAX_DISTANCE_FROM_LINEAR
            :dedent: 4

    .. md-tab-item:: Python

        .. literalinclude:: ../../../python/tests/test_ccd.py
            :language: python
            :start-after: # BEGIN_RIGID_2D_MAX_DISTANCE_FROM_LINEAR
            :end-before: # END_RIGID_2D_MAX_DISTANCE_FROM_LINEAR
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

Interval-Based Nonlinear CCD
----------------------------

.. warning::
    The following section requires enabling the ``filib`` interval arithmetic library. ``filib`` is licensed under L-GPL 2.1, so special care must be taken when using it. See the `filib dependency note <../../cpp.html#filib-dependency-note>`_ for more information.

If an analytic expression for the ``max_distance_from_linear`` function is not available, we can use interval arithmetic to compute a conservative envelope.

.. md-tab-set::

    .. md-tab-item:: C++

        The :cpp:class:`ipc::IntervalNonlinearTrajectory` class does this for us and all we need to provide is a function to compute the point's position over a time interval.

    .. md-tab-item:: Python

        The :py:class:`ipctk.IntervalNonlinearTrajectory` class does this for us and all we need to provide is a function to compute the point's position over a time interval.

Conservative Envelope
~~~~~~~~~~~~~~~~~~~~~

Our implementation of the ``max_distance_from_linear`` function is as follows.

Start by defining a linear interpolation function:

.. math::
    \operatorname{lerp}(a, b, t) := (b - a) t + a,

which interpolates between two points :math:`a` and :math:`b` at time :math:`t`.

The exact envelope from above is bounded by a interval approximation:

.. math::
    \begin{align}
    &\max_{t \in [0, 1]} \| p(\operatorname{lerp}(t_0, t_1, t)) - \operatorname{lerp}(p(t_0), p(t_1), t) \|_2\\
    &\leq \sup(\| p_{\Box}([t_0, t_1]) - ((p(t_1) - p(t_0)) \cdot [0, 1] + p(t_0)) \|_2).
    \end{align}

Therefore, we can compute the a conservative estimate of the envelope by implementing :math:`p_{\Box}([t_0, t_1])` and :math:`p(t)`.

.. note::
    In practice we subdivide the interval into smaller intervals and compute the envelope over each subinterval. This is done to create a more accurate estimate.

Interval Arithmetic
~~~~~~~~~~~~~~~~~~~

`Interval arithmetic <https://en.wikipedia.org/wiki/Interval_arithmetic>`_ is a method to compute bounds on the range of a function over an interval. We can construct a vector of intervals to represent the position of the point over a time interval.

The following code snippet shows an example of how to use interval arithmetic to compute the position of a point over a time interval:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: cpp

            #include <ipc/utils/interval.hpp>

            using namespace ipc;

            Vector2I position(
                const VectorMax3d& center,
                const VectorMax3d& point,
                const double omega,
                const Interval& t)
            {
                // 2×2 matrix of intervals representing the rotation matrix
                Matrix2I R;
                R << cos(omega * t), -sin(omega * t),
                     sin(omega * t),  cos(omega * t);
                return R * (point - center) + center;
            }

        The full documentation for the ``Interval`` class can be found `in the C++ API <../../cpp-api/interval.html>`_.

    .. md-tab-item:: Python

        .. code-block:: python

            import numpy as np

            from ipctk.filib import Interval, sin, cos

            def position(center, point, omega, t : Interval):
                R = np.array([
                    [cos(omega * t), -sin(omega * t)],
                    [sin(omega * t),  cos(omega * t)]
                ])
                return R @ (point - center) + center

        The full documentation for the filib python bindings can be found `in the Python API <../python-api/interval.html>`_.
