Advanced Friction
=================

This tutorial covers some advanced features of friction in the IPC Toolkit, including
spatially varying coefficients of friction and separate coefficients of friction for static and kinetic (dynamic) friction.

.. seealso::

    For an introduction to friction modeling and the dissipative potential, see
    `"Getting Started" <getting_started.html#friction>`_.


Spatially Varying Coefficients of Friction
------------------------------------------

.. warning::
    :title: Deprecated

    Spatially varying coefficient of friction is achieved by assigning coefficients to each vertex in the mesh. However, friction coefficients are not a material property and should instead be assigned to the contact pair. This feature will be replaced with a per-pair friction coefficient in a future release.

You can specify spatially varying coefficients of friction by passing an ``Eigen::VectorXd`` to ``TangentialCollisions::build``. Each entry in the vector corresponds to the coefficient of friction for a specific vertex in the mesh. This allows you to assign different friction coefficients to different parts of the mesh, enabling more realistic simulations of complex materials and surfaces.

You can also provide an optional ``blend_mu`` parameter to blend the coefficient of friction on either side of the contact. The default behavior is to average the coefficients of friction on both sides, but you can specify a custom blending function if needed (e.g., multiplying them or taking the maximum or minimum).

Separate Coefficients for Static and Kinetic Friction
-----------------------------------------------------

.. tip::
    :title: Experimental Feature

    The support for separate static and kinetic friction coefficients is an experimental feature introduced in v1.4.0. The behavior may be improved in future releases.

Traditionally, friction is modeled with a single coefficient. However, physical systems often exhibit different resistance to motion when at rest (static friction) versus when sliding (kinetic friction). The IPC Toolkit supports specifying both:

- **Static friction coefficient** (``mu_s``): Governs the maximum force resisting the initiation of motion.
- **Kinetic friction coefficient** (``mu_k``): Governs the force resisting ongoing motion.

Usage
~~~~~

To use separate static and kinetic friction coefficients, you can pass them as parameters when building the tangential collisions:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            ipc::TangentialCollisions tangential_collisions;
            tangential_collisions.build(
                collision_mesh, vertices, collisions, B, barrier_stiffness,
                mu_s, mu_k);

    .. md-tab-item:: Python

        .. code-block:: python

            tangential_collisions = ipctk.TangentialCollisions()
            tangential_collisions.build(
                collision_mesh, vertices, collisions, B, barrier_stiffness,
                mu_s, mu_k)

How It Works
~~~~~~~~~~~~

Background
^^^^^^^^^^

IPC defines the local friction force as

.. math::
   F = -\mu \lambda \mathbf{T}(x) f_1(\|\mathbf{u}\|) \frac{\mathbf{u}}{\|\mathbf{u}\|},

where :math:`\lambda` is the contact force magnitude, :math:`T(x) \in \mathbb{R}^{3 n \times 2}` is the oriented sliding basis, and :math:`\mu` is the local friction coefficient. The function :math:`f_1` is a mollifier that smooths the transition between static and kinetic friction defined as

.. math::
   f_1(y)= \begin{cases}
        -\frac{y^2}{\epsilon_v^2}+\frac{2 y}{\epsilon_v} & \text{for}\: 0 \leq y \leq \epsilon_v \\
        1 & \text{otherwise}
    \end{cases}

where :math:`\epsilon_v` is a small constant (e.g., ``0.001``). The following plot show the behavior of the function :math:`f_1`:

.. figure:: ../_static/img/f1.png
   :align: center

   Smooth mollifier :math:`f_1(\|\mathbf{u}\|)` where :math:`\epsilon_v = 0.001`.

To create a dissipative potential we integrate :math:`f_1` to obtain a smooth mollifier :math:`f_0`:

.. math::
    f_0(y)= \begin{cases}
        \frac{\epsilon_{v}}{3} + \frac{y^{2}}{\epsilon_{v}} - \frac{y^{3}}{3 \epsilon_{v}^{2}} & \text{for}\: 0 \leq y \leq \epsilon_{v} \\
        y & \text{otherwise} \\
    \end{cases}


The following plot show the behavior of the function :math:`f_0`:

.. figure:: /_static/img/f0.png
   :align: center

   Integrated mollifier :math:`f_0(\|\mathbf{u}\|)` where :math:`\epsilon_v = 0.001`.

Smooth :math:`\mu`
^^^^^^^^^^^^^^^^^^

When adding separate coefficients for static and kinetic friction, we need to maintain the :math:`C^1` continuity of the friction force. This lead us to define a smooth coefficient of friction :math:`\mu(y)` that transitions between the static and kinetic coefficients based on the magnitude of the relative velocity :math:`y = \|\mathbf{u}\|`. The smooth coefficient of friction is defined as

.. math::
    \mu(y) = \begin{cases} \mu_{s} + \left(\mu_{k} - \mu_{s}\right) \left(\frac{3 y^{2}}{\epsilon_{v}^{2}} - \frac{2 y^{3}}{\epsilon_{v}^{3}}\right) & \text{for}\: 0 \leq y \leq \epsilon_{v}\\
        \mu_{k} & \text{otherwise}.
    \end{cases}

We plot the smooth coefficient of friction :math:`\mu(y)` below:

.. figure:: /_static/img/mu.png
   :align: center

   Smooth coefficient of friction :math:`\mu(\|\mathbf{u}\|)` where :math:`\mu_s = 1`, :math:`\mu_k = 0.1`, and :math:`\epsilon_v = 0.001`.

Smooth :math:`\mu` Mollifier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Replacing the constant coefficient of friction :math:`\mu` with a smooth function :math:`\mu(\|\mathbf{u}\|)` allows us to smoothly transition between static and kinetic friction. The function :math:`\mu(\|\mathbf{u}\|) f_1(\|\mathbf{u}\|)` is plotted below:

.. figure:: /_static/img/mu_f1.png
   :align: center

   Smooth coefficient of friction multiplied by the friction mollifier :math:`\mu(\|\mathbf{u}\|) f_1(\|\mathbf{u}\|)` where :math:`\mu_s = 1`, :math:`\mu_k = 0.1`, and :math:`\epsilon_v = 0.001`.

However, :math:`\frac{\mathrm{d}}{\mathrm{d}y} \mu(y) f_0(y) \neq \mu(y) f_1(y)`, so we need to adjust the integrated mollifier by integrating the product of the smooth coefficient of friction and the mollifier:

.. math::
    \int \mu(y) f_1(y) \mathrm{d} y = \begin{cases}
        \frac{\epsilon_{v}^{6} \left(17 \mu_{k} - 7 \mu_{s}\right) + 30 \epsilon_{v}^{4} \mu_{s} y^{2} - 10 \epsilon_{v}^{3} \mu_{s} y^{3} + 45 \epsilon_{v}^{2} \left(\mu_{k} - \mu_{s}\right) y^{4} + 42 \epsilon_{v} \left(- \mu_{k} + \mu_{s}\right) y^{5} + 10 \left(\mu_{k} - \mu_{s}\right) y^{6}}{30 \epsilon_{v}^{5}} & \text{for}\: \epsilon_{v} > y \\
        \mu_{k} y & \text{otherwise}
    \end{cases}.




The following plot shows the behavior of the integrated mollifier multiplied by the smooth coefficient of friction:

.. figure:: /_static/img/int_mu_f1_dx.png
   :align: center

   Integrated mollifier :math:`\int \mu(y) f_1(y) \mathrm{d} y` where :math:`\mu_s = 1`, :math:`\mu_k = 0.1`, and :math:`\epsilon_v = 0.001`.

Discussion
^^^^^^^^^^

While this approach provides a smooth transition between static and kinetic friction, there are some considerations:

1. The product :math:`\mu(y) f_1(y)` underestimates the friction force in the static regime, which may lead to less accurate simulations of static friction.
    - This is, when :math:`\mu_k < \mu_s` and :math:`|y| \leq \epsilon_v`, :math:`\max(\mu(y) f_1(y)) < \mu_s`.
    - We could address this by scaling by :math:`\frac{\max(\mu(y) f_1(y))}{\mu_s}`, but computing the maximum is non-trivial.
2. The combined function :math:`\mu(y) f_1(y)` is a degree 5 polynomial, which is more complex than the original mollifier :math:`f_1(y)` (degree 2). This may lead to more difficult to optimize problems. There are two options to address this:
    a. Replacing :math:`\mu(y)` with a piecewise quadratic function would reduce the degree to 4, but it would still be more complex than the original mollifier.
    b. Replacing :math:`\mu(y) f_1(y)` with a piecewise quadratic function that has the desired behavior. This help with 1. as well. However making it backwards compatible with the original mollifier is challenging.

If you have suggestions for improving this approach or alternative methods, please reach out on our `GitHub Discussions <https://github.com/ipc-sim/ipc-toolkit/discussions>`.

Future Directions
-----------------

The IPC Toolkit is continuously evolving, and future releases may include:

- Anisotropic friction models that account for direction-dependent friction.
- Velocity-dependent friction models that adjust friction coefficients based on relative velocity magnitude.
- Rolling coefficients of friction for scenarios involving rolling contacts.

We encourage community contributions to expand these advanced friction models. Feel free to submit pull requests with your improvements or open a discussion on GitHub to propose new features.