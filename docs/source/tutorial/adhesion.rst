Adhesion
========

Adhesion simulates the sticking interaction between surfaces, such as glue bonding or material contact.

We provide functionality to compute normal adhesion (perpendicular forces) and tangential adhesion (parallel forces) using the model of :cite:t:`Fang2023AugmentedStickyInteractions`.

Normal Adhesion
---------------

The normal adhesion potential models the attraction force based on the distance between two surfaces.

.. math::
    A_n(\mathbf{x}) = \sum_{k \in C} a(d_k(\mathbf{x}); \hat{d}_p, \hat{d}_a, Y, \epsilon_c)

For a given distance :math:`d`:

.. math::
    a(d)= \begin{cases}
    a_2\left(1-\frac{\hat{d}_a}{\hat{d}_p}\right) d^2+a_2\left(\hat{d}_a-\hat{d}_p\right)^2-d_p^2 a_1 & 0<d<\hat{d}_p \\
    a_2 d^2 - 2 a_2 \hat{d}_a d+a_2 \hat{d}_a^2 & \hat{d}_p \leq d<\hat{d}_a \\
    0 & \hat{d}_a \leq d
    \end{cases}

and

.. math::
    a_2=\frac{Y \varepsilon_c}{2\left(\hat{d}_p-\hat{d}_a\right)}.

The parameters are:

- :math:`\hat{d}_p` (``dhat_p``) is the threshold distance for adhesion (in units of distance) where the largest adhesion force is applied
- :math:`\hat{d}_a` (``dhat_a``) is the adhesion activation distance (in units of distance), representing the maximum range where adhesion forces are active
- :math:`Y` (``Y``) is the adhesion stiffness (in units of stress, such as Young's modulus), controlling the intensity of adhesion forces
- :math:`\epsilon_c` (``eps_c``) is the adhesion coefficient (unitless) that defines the critical strain at which adhesion forces decrease

We can build a normal adhesion potential object and compute the adhesion potential for a given set of normal collisions.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

                const double dhat_p = 1e-3;
                const double dhat_a = 2 * dhat_p;
                const double Y = 1e3;
                const double eps_c = 0.5;

                const ipc::NormalAdhesionPotential A_n(dhat_p, dhat_a, Y, eps_c)
                double adhesion_potential = A_n(normal_collisions, collision_mesh, vertices);

    .. md-tab-item:: Python

        .. code-block:: python

                dhat_p = 1e-3
                dhat_a = 2 * dhat_p
                Y = 1e3
                eps_c = 0.5

                A_n = ipctk.NormalAdhesionPotential(dhat_p, dhat_a, Y, eps_c)
                adhesion_potential = A_n(normal_collisions, collision_mesh, vertices)

Derivatives
^^^^^^^^^^^

We can also compute the first and second derivatives of the normal adhesion potential with respect to the vertex positions.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::VectorXd adhesion_potential_grad =
                A_n.gradient(normal_collisions, collision_mesh, vertices);

            Eigen::SparseMatrix<double> adhesion_potential_hess =
                A_n.hessian(normal_collisions, collision_mesh, vertices);

    .. md-tab-item:: Python

        .. code-block:: python

            adhesion_potential_grad = A_n.gradient(normal_collisions, collision_mesh, vertices)

            adhesion_potential_hess = A_n.hessian(normal_collisions, collision_mesh, vertices)

Tangential Adhesion
-------------------

The tangential adhesion potential models resistance to sliding (parallel to surfaces).

It is structured similar to the friction model with a smooth transition between sticking and sliding and lagged normal forces and tangential bases.

.. math::
    A_{t}(\mathbf{u}) = \sum_{k \in C} \mu_a \lambda_{k}^a f_{0}^a(\|\mathbf{T}_k^⊤ \mathbf{u}\|; \epsilon_a)

where :math:`C` is the lagged collisions, :math:`\lambda_k^a` is the normal adhesion force magnitude for the :math:`k`-th collision, :math:`\mathbf{T}_k` is the tangential basis for the :math:`k`-th collision, and :math:`f_0^a` is the smooth tangential adhesion function used to approximate the non-smooth transition from sticking to sliding.

For relative displacement magnitude :math:`y`:

.. math::
    f_0^a(y) = \begin{cases}
    -\frac{y^3}{3\epsilon_a^2} + \frac{y^2}{\epsilon_a} & 0 < y < 2 \epsilon_a \\
    \frac{4 \epsilon_a}{3} & 2 \epsilon_a \leq y
    \end{cases}

where :math:`\epsilon_a` (``eps_a``) is the adhesion threshold (in units of displacement) used to smoothly transition from sticking to sliding.

We can build a tangential adhesion potential object and compute the adhesion potential for a given set of tangential collisions.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            const double eps_a = 0.01;
            const ipc::TangentialAdhesionPotential A_t(eps_a);
            double adhesion_potential = A_t(tangential_collisions, collision_mesh, displacement);

    .. md-tab-item:: Python

        .. code-block:: python

            eps_a = 0.01
            A_t = ipctk.TangentialAdhesionPotential(eps_a)
            adhesion_potential = A_t(tangential_collisions, collision_mesh, displacement);

Derivatives
^^^^^^^^^^^

We can also compute the first and second derivatives of the tangential adhesion potential with respect to the displacement.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::VectorXd adhesion_potential_grad =
                A_t.gradient(tangential_collisions, collision_mesh, displacement);

            Eigen::SparseMatrix<double> adhesion_potential_hess =
                A_t.hessian(tangential_collisions, collision_mesh, displacement);

    .. md-tab-item:: Python

        .. code-block:: python

            adhesion_potential_grad = A_t.gradient(tangential_collisions, collision_mesh, displacement)

            adhesion_potential_hess = A_t.hessian(tangential_collisions, collision_mesh, displacement)
