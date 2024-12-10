Adhesion
--------
Adhesion simulates the sticking interaction between surfaces, such as glue bonding or material contact. 
Provides functionality to compute normal adhesion (perpendicular forces) and tangential adhesion (parallel forces).

Normal Adhesion
^^^^^^^^^^^^^^^

For a given distance :math:`d`:

- If :math:`d < \hat{d}_p`:
  
  .. math::
     A_n(d) = a_1 \cdot d^2 + c_1, \quad \text{where } a_1 = a_2 \left(1 - \frac{\hat{d}_a}{\hat{d}_p}\right)

- If :math:`\hat{d}_p \leq d < \hat{d}_a`:

  .. math::
     A_n(d) = (a_2 \cdot d + b_2) \cdot d + c_2, \quad \text{where } b_2 = -2 a_2 \hat{d}_a, \, c_2 = a_2 \hat{d}_a^2

- If :math:`d \geq \hat{d}_a`:

  .. math::
     A_n(d) = 0

The normal adhesion potential models the attraction force based on the distance between two surfaces.

Here dhat_p (:math:\hat{d}_p) is the threshold distance for adhesion (in units of distance) where the largest adhesion force is applied.
Here dhat_a (:math:\hat{d}_a) is the adhesion activation distance (in units of distance), representing the maximum range where adhesion forces are active.
Here Y (:math:Y) is the adhesion stiffness (in units of stress, such as Young's modulus), controlling the intensity of adhesion forces.
Here eps_c (:math:\epsilon_c) is the adhesion coefficient (unitless, critical strain) that defines the critical strain at which adhesion forces decrease.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

                const double dhat_p = 1e-3;
                const double dhat_a = 2 * dhat_p;
                const double Y = 1e3;
                const double eps_c = 0.5;

                const NormalAdhesionPotential A(dhat_p, dhat_a, Y, eps_c)
                double adhesion_potential = A(normal_collisions, collision_mesh, velocity);

    .. md-tab-item:: Python

        .. code-block:: python

                dhat_p = 1e-3
                dhat_a = 2 * dhat_p
                Y = 1e3
                eps_c = 0.5

                A = NormalAdhesionPotential(dhat_p, dhat_a, Y, eps_c)
                adhesion_potential = A(normal_collisions, collision_mesh, velocity)

Normal Derivatives
^^^^^^^^^^^

We can also compute the first and second derivatives of the normal adhesion potential with respect to the velocities.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::VectorXd adhesion_potential_grad =
                A.gradient(normal_collisions, collision_mesh, velocity);

            Eigen::SparseMatrix<double> adhesion_potential_hess =
                A.hessian(normal_collisions, collision_mesh, velocity);

    .. md-tab-item:: Python

        .. code-block:: python

            adhesion_potential_grad = A.gradient(
                normal_collisions, collision_mesh, velocity)

            adhesion_potential_hess = A.hessian(
                normal_collisions, collision_mesh, velocity)

Tangential Adhesion
^^^^^^^^^^^^^^^

The tangential adhesion potential models resistance to sliding (parallel to surfaces).

For velocity :math:`y`:

- If :math:`0 \leq y < 2 \varepsilon_a`:

  .. math::
     A_t(y) = \frac{y^2}{\varepsilon_a} \left(1 - \frac{y}{3 \varepsilon_a}\right)

- If :math:`y \geq 2 \varepsilon_a`:

  .. math::
     A_t(y) = \frac{4 \varepsilon_a}{3}

Here ``eps_a`` (:math:`\epsilon_a`) is the adhesion threshold (in units of velocity) used to smoothly transition.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            const double eps_a = 0.01;
            const TangentialAdhesionPotential A(eps_a);
            double adhesion_potential = A(tangential_collisions, collision_mesh, velocity);
    
    .. md-tab-item:: Python

        .. code-block:: python

            eps_a = 0.01
            A = TangentialAdhesionPotential(eps_a)
            adhesion_potential = A(tangential_collisions, collision_mesh, velocity);

Derivatives
^^^^^^^^^^^

We can also compute the first and second derivatives of the tangential adhesion potential with respect to the velocities.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::VectorXd adhesion_potential_grad =
                A.gradient(tangential_collisions, collision_mesh, velocity);

            Eigen::SparseMatrix<double> adhesion_potential_hess =
                A.hessian(tangential_collisions, collision_mesh, velocity);

    .. md-tab-item:: Python

        .. code-block:: python

            adhesion_potential_grad = A.gradient(
                tangential_collisions, collision_mesh, velocity)

            adhesion_potential_hess = A.hessian(
                tangential_collisions, collision_mesh, velocity)
