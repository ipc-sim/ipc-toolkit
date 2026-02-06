.. _geometric-contact-potential-formulation:

Geometric Contact Potential
===========================

In addition to the original implementation of :cite:t:`Li2020IPC`, we also implement the Geometric Contact Potential (GCP) from :cite:t:`Huang2025GCP`.

GCP is implemented as separate collision and potential classes. A basic example of computing the GCP potential is as follows.

.. code-block:: c++

    ipc::SmoothCollisions collisions;
    if (use_adaptive_dhat)
        collisions.compute_adaptive_dhat(collision_mesh, vertices, params);
    collisions.build(collision_mesh, vertices, params, use_adaptive_dhat);

    ipc::SmoothContactPotential barrier_potential;
    double b = barrier_potential(collisions, mesh, vertices);

.. important::
    If `use_adaptive_dhat` is true, make sure to call `SmoothCollisions::compute_adaptive_dhat()` before `SmoothCollisions::build()`.

Technical Details
-----------------

Similar to IPC, GCP potential has a distance barrier that converges to infinity as the distance converges to zero.

.. math::
    B(d; \hat{d}, p) := \frac{3}{2} B^3(\frac{2d}{\hat{d}}) d^{-p}

where :math:`\hat{d}` is the support size of the barrier function, :math:`p=d-1`, where :math:`d` is the dimension of the physics scene (2 or 3).

The main motivation of developing GCP is to eliminate the spurious collision forces in IPC when :math:`\hat{d}` is relatively large. The two important components to achieve this are *local minimum constraint* and *exterior direction constraint*. Please check out :cite:t:`Huang2025GCP` for more details, here we only give a high level idea.

The local minimum constraint is by noticing that the barrier potential is only necessary when the two primitives in a collision pair are the "distance local minima" on each surface they live. To achieve this, we consider the gradient of distance function along the tangent directions on the surface, and smoothly clamp the barrier to zero if the gradient of distance is positive along the tangent directions, since it suggests that there already exists a closer point on the surface that can trigger the barrier.

The exterior direction constraint is by noticing that, for a thin volumetric shell, the barrier potential is not needed for points on the opposite sides of the shell. To achieve this, we smoothly clamp the barrier based on the normal directions of two primitives in a collision pair. If the normal direction points to the opposite direction of the distance direction, we can skip this collision pair.

Both constraints not only reduces the number of collision pairs, but also increases the accuracy of collision handling by eliminating spurious forces that are non-physical, allowing wider range of :math:`\hat{d}` to be used in the simulation.

Parameter choices
-----------------

Besides :math:`\hat{d}`, GCP also has other parameters to control the smoothness of local minimum and exterior direction constraints. For each constraint, there are :math:`\alpha` and :math:`\beta` to control the support region of the smooth Heaviside function (:math:`\alpha_t,\ \beta_t` for local minimum constraint, and :math:`\alpha_n,\ \beta_n` for exterior direction constraint). The smooth Heaviside function :math:`H(z;\alpha,\beta)\in C^1(\mathbb{R})` satisfies that

.. math::
    H(z;\alpha,\beta) = 0,\ \forall z < -\alpha \\
    H(z;\alpha,\beta) = 1,\ \forall z > \beta \\
    H'(z;\alpha,\beta) \geq 0,\ \forall z

Since :math:`H(z;\alpha,\beta)` takes the dot and cross products of unit vectors as inputs, we have :math:`-1\leq z\leq 1`. Thus, the basic requirement for parameters is:

.. math::
    0 \leq -\beta < \alpha \leq 1.

As :math:`\alpha + \beta` decreases, the support size of the Heaviside function reduces and the function becomes less smooth, thus, the overall potential becomes harder to optimize.

.. tip::
    For simplicity, we recommend to just use :math:`\beta_n = \beta_t = 0`. To make sure the potential landscape is smooth, the recommended parameter choices are :math:`\alpha_t\in [0.2, 0.9]`, :math:`\alpha_n = 0.1`.

Friction
--------

We implement the original friction formulation of :cite:t:`Li2020IPC` following the same style, details are covered in the paper :cite:t:`Huang2025GCP`.
