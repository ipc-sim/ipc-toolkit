Getting Started
===============

This tutorial will walk you through the basics of using the IPC Toolkit.

Collision Mesh
--------------

First, we need to create a collision mesh. The ``CollisionMesh`` data structure represents the surface geometry used for collision processing.

We will start by creating a collision mesh from a ``bunny.obj`` mesh file (you can find the mesh `here <https://github.com/ipc-sim/ipc-toolkit/blob/main/tests/data/bunny.obj>`_):

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <ipc/ipc.hpp>

            #include <Eigen/Core>
            #include <Eigen/Sparse>
            #include <igl/readOBJ.h>
            #include <igl/edges.h>

            // ...

            Eigen::MatrixXd rest_positions;
            Eigen::MatrixXi edges, faces;
            igl::readOBJ("bunny.obj", rest_positions, faces);
            igl::edges(faces, edges);

            ipc::CollisionMesh collision_mesh(rest_positions, edges, faces);

    .. md-tab-item:: Python

        .. code-block:: python

            import ipctk

            import meshio

            mesh = meshio.read("bunny.obj")
            rest_positions = mesh.points
            faces = mesh.cells[0].data
            edges = ipctk.edges(faces)

            collision_mesh = ipctk.CollisionMesh(rest_positions, edges, faces)

The matrix ``rest_positions`` contains the undeformed positions of the vertices. It should be :math:`|V| \times d` where :math:`|V|` is the number of vertices and :math:`d \in \{2, 3\}` is the dimension.
Each row of the ``edges`` and ``faces`` matrices contain the vertex IDs (row number in ``rest_positions``) of the edge or face's vertices.
The size of ``edges`` and ``faces`` are ``#E x 2`` and ``#F x 3`` respectively (``#E`` and ``#F`` are the number of edges and faces).

.. note::
   Only linear triangular faces are supported. If your mesh has nonlinear or non-triangular faces, you will need to triangulate them.

.. note::
   In 2D only the ``edges`` matrix is required. In 3D both ``edges`` and ``faces`` are required.

Collision Constraints
---------------------

Now that we have a collision mesh, we can compute the contact potential. To do this we first need to build the set of active collision constraints (``CollisionConstraints``).

To start we need the current positions of the ``vertices``. For this tutorial, let us use squash the bunny has to 1% of its original height.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::MatrixXd vertices = collision_mesh.rest_positions();
            vertices.col(1) *= 0.01; // Squash the bunny in the y-direction

    .. md-tab-item:: Python

        .. code-block:: python

            vertices = collision_mesh.rest_positions()
            vertices[:, 1] *= 0.01  # Squash the bunny in the y-direction

Using these deformed positions, we can build the set of active collision constraints.
For this, we need to specify :math:`\hat{d}` (``dhat``), the IPC barrier activation distance.
We will use a value of :math:`\hat{d} = 10^{-3}`.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            const double dhat = 1e-3;

            ipc::CollisionConstraints collision_constraints;
            collision_constraints.build(collision_mesh, vertices, dhat);

    .. md-tab-item:: Python

        .. code-block:: python

            dhat = 1e-3

            collision_constraints = ipctk.CollisionConstraints()
            collision_constraints.build(collision_mesh, vertices, dhat)

This will automatically use a spatial data structure to perform a broad-phase culling and then perform a narrow-phase culling by computing distances (discarding any constraint with a distance :math:`> \hat{d}`).

Contact Potential
^^^^^^^^^^^^^^^^^

Now we can compute the contact potential using the ``CollisionConstraints``.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            double contact_potential = collision_constraints.compute_potential(
                collision_mesh, vertices, dhat);

    .. md-tab-item:: Python

        .. code-block:: python

            contact_potential = collision_constraints.compute_potential(
                collision_mesh, vertices, dhat)

This returns a scalar value ``contact_potential`` which is the sum of the contact potential for each active constraint.

Mathematically this is defined as

.. math::
   B(x) = \sum_{k \in C} b\left(d_k(x)\right),

where :math:`C` is the active collision constraints, :math:`d_k` is the distance (squared) of the :math:`k`-th active constraint, and :math:`b` is IPC's C2-clamped log-barrier function.

.. note::
   This is not premultiplied by the barrier stiffness :math:`\kappa`.

Contact Potential Derivative
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We can also compute the first and second derivatives of the contact potential with respect to the vertex positions.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::VectorXd contact_potential_grad =
                collision_constraints.compute_potential_gradient(
                    collision_mesh, vertices, dhat);

            Eigen::SparseMatrix<double> contact_potential_hess =
                collision_constraints.compute_potential_hessian(
                    collision_mesh, vertices, dhat);

    .. md-tab-item:: Python

        .. code-block:: python

            contact_potential_grad = collision_constraints.compute_potential_gradient(
                collision_mesh, vertices, dhat)

            contact_potential_hess = collision_constraints.compute_potential_hessian(
                collision_mesh, vertices, dhat)

These return the gradient and hessian of the contact potential as a dense vector and sparse matrix respectively.

The derivatives are taken with respect to the row-wise flattened vertices. That is, for ``vertices``

.. math::
    \begin{bmatrix}
    x_1 & y_1 & z_1 \\
    & \vdots & \\
    x_n & y_n & z_n \\
    \end{bmatrix},

you will get the gradient of size :math:`|V|d \times 1` with the order

.. math::
    \begin{bmatrix}
    \frac{\partial B}{\partial x_1} &
    \frac{\partial B}{\partial y_1} &
    \frac{\partial B}{\partial z_1} &
    \cdots &
    \frac{\partial B}{\partial x_n} &
    \frac{\partial B}{\partial y_n} &
    \frac{\partial B}{\partial z_n}
    \end{bmatrix}^T,

and the Hessian of size :math:`|V|d \times |V|d` with the order

.. math::
    \begin{bmatrix}
    \frac{\partial^2 B}{\partial x_1^2} &
    \frac{\partial^2 B}{\partial x_1 \partial y_1} &
    \frac{\partial^2 B}{\partial x_1 \partial z_1} &
    \cdots &
    \frac{\partial^2 B}{\partial x_1 \partial x_n} &
    \frac{\partial^2 B}{\partial x_1 \partial y_n} &
    \frac{\partial^2 B}{\partial x_1 \partial z_n} \\
    %
    \frac{\partial^2 B}{\partial y_1 \partial x_1} &
    \frac{\partial^2 B}{\partial y_1^2} &
    \frac{\partial^2 B}{\partial y_1 \partial z_1} &
    \cdots &
    \frac{\partial^2 B}{\partial y_1 \partial x_n} &
    \frac{\partial^2 B}{\partial y_1 \partial y_n} &
    \frac{\partial^2 B}{\partial y_1 \partial z_n} \\
    %
    \frac{\partial^2 B}{\partial z_1 \partial x_1} &
    \frac{\partial^2 B}{\partial z_1 \partial y_1} &
    \frac{\partial^2 B}{\partial z_1^2} &
    \cdots &
    \frac{\partial^2 B}{\partial z_1 \partial x_n} &
    \frac{\partial^2 B}{\partial z_1 \partial y_n} &
    \frac{\partial^2 B}{\partial z_1 \partial z_n} \\
    %
    \vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
    %
    \frac{\partial^2 B}{\partial x_n \partial x_1} &
    \frac{\partial^2 B}{\partial x_n \partial y_1} &
    \frac{\partial^2 B}{\partial x_n \partial z_1} &
    \cdots &
    \frac{\partial^2 B}{\partial x_n^2} &
    \frac{\partial^2 B}{\partial x_n \partial y_n} &
    \frac{\partial^2 B}{\partial x_n \partial z_n} \\
    %
    \frac{\partial^2 B}{\partial y_n \partial x_1} &
    \frac{\partial^2 B}{\partial y_n \partial y_1} &
    \frac{\partial^2 B}{\partial y_n \partial z_1} &
    \cdots &
    \frac{\partial^2 B}{\partial y_n \partial x_n} &
    \frac{\partial^2 B}{\partial y_n^2} &
    \frac{\partial^2 B}{\partial y_n \partial z_n} \\
    %
    \frac{\partial^2 B}{\partial z_n \partial x_1} &
    \frac{\partial^2 B}{\partial z_n \partial y_1} &
    \frac{\partial^2 B}{\partial z_n \partial z_1} &
    \cdots
    &
    \frac{\partial^2 B}{\partial z_n \partial x_n} &
    \frac{\partial^2 B}{\partial z_n \partial y_n} &
    \frac{\partial^2 B}{\partial z_n^2}
    \end{bmatrix}.

Adaptive Barrier Stiffness
^^^^^^^^^^^^^^^^^^^^^^^^^^

The last piece of the contact potential is the barrier stiffness. This is a weight that is multiplied by the barrier potential to better scale it relative to the energy potential. This can be a fixed value or adaptive.

To compute the adaptive barrier stiffness, we can use two functions: ``initial_barrier_stiffness`` and ``update_barrier_stiffness``. The function ``initial_barrier_stiffness`` compute the initial value from the current energy and contact potential gradients. This function also provides a minimum and maximum value for the barrier stiffness. The function ``update_barrier_stiffness`` updates the barrier stiffness if the minimum distance has become too small.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // (beginning of nonlinear solve)

            Eigen::VectorXd grad_energy = ...; // gradient of elastic energy potential
            Eigen::VectorXd grad_contact = collision_constraints.compute_potential_gradient(
                collision_mesh, vertices, dhat);

            double bbox_diagonal = ipc::world_bbox_diagonal_length(vertices);

            double max_barrier_stiffness; // output of initial_barrier_stiffness
            double barrier_stiffness = ipc::initial_barrier_stiffness(
                bbox_diagonal, dhat, avg_mass, grad_energy, grad_contact,
                max_barrier_stiffness);

            double prev_distance =
                collision_constraints.compute_minimum_distance(
                    collision_mesh, vertices);

            // ...

            // (end of nonlinear iteration)

            double curr_distance =
                collision_constraints.compute_minimum_distance(collision_mesh, vertices);

            barrier_stiffness = ipc::update_barrier_stiffness(
                prev_distance, curr_distance, max_barrier_stiffness, barrier_stiffness,
                bbox_diagonal);

            prev_distance = curr_distance;

            // (next iteration)

    .. md-tab-item:: Python

        .. code-block:: python

            # (beginning of nonlinear solve)

            grad_energy = ...  # gradient of elastic energy potential
            grad_contact = collision_constraints.compute_potential_gradient(
                collision_mesh, vertices, dhat)

            bbox_diagonal = ipctk.world_bbox_diagonal_length(vertices)

            barrier_stiffness, max_barrier_stiffness = ipctk.initial_barrier_stiffness(
                bbox_diagonal, dhat, avg_mass, grad_energy, grad_contact,
                max_barrier_stiffness)

            prev_distance = collision_constraints.compute_minimum_distance(
                collision_mesh, vertices)

            # ...

            # (end of nonlinear iteration)

            curr_distance = collision_constraints.compute_minimum_distance(
                collision_mesh, vertices)

            barrier_stiffness = ipctk.update_barrier_stiffness(
                prev_distance, curr_distance, max_barrier_stiffness, barrier_stiffness,
                bbox_diagonal)

            prev_distance = curr_distance

            # (next iteration)

Friction
--------

Computing the friction dissipative potential is similar to the contact potential, but because it is a lagged model, we need to build it from a fixed set of constraints.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            ipc::FrictionConstraints friction_constraints;
            collision_constraints.build(
                collision_mesh, vertices, contact_constraints, dhat, barrier_stiffness, mu);

    .. md-tab-item:: Python

        .. code-block:: python

            friction_constraints = ipctk.FrictionConstraints()
            friction_constraints.build(
                collision_mesh, vertices, contact_constraints, dhat, barrier_stiffness, mu)

Here ``mu`` (:math:`\mu`) is the (global) coefficient of friction, and ``barrier_stiffness`` (:math:`\kappa`) is the barrier stiffness.

Friction Dissipative Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now we can compute the friction dissipative potential using the ``FrictionConstraints``.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            double friction_potential = friction_constraints.compute_potential(
                collision_mesh, velocity, epsv);

    .. md-tab-item:: Python

        .. code-block:: python

            friction_potential = friction_constraints.compute_potential(
                collision_mesh, velocity, epsv)

Here ``epsv`` (:math:`\epsilon_v`) is the static friction threshold (in units of velocity) used to smoothly transition from dynamic to static friction.

.. important::
   The friction potential is a function of the velocities rather than the positions. We can compute the velocities directly from the current and previous position(s) based on our time-integration scheme. For example, if we are using backward Euler integration, then the velocity is

   .. math::
      v = \frac{x - x_{n-1}}{h},

   where :math:`x` is the current position, :math:`x_{n-1}` is the previous position, and :math:`h` is the time step size.

This returns a scalar value ``friction_potential`` which is the sum of the individual friction potentials.

Mathematically this is defined as

.. math::
   D(x) = \sum_{k \in C} \mu\lambda_k^nf_0\left(\|T_k^Tv\|; \epsilon_v\right),

where :math:`C` is the lagged collision constraints, :math:`\lambda_k^n` is the normal force magnitude for the :math:`k`-th contact, :math:`T_k` is the tangential basis for the :math:`k`-th contact, and :math:`f_0` is the smooth friction function used to approximate the non-smooth transition from dynamic to static friction.

Derivatives
^^^^^^^^^^^

We can also compute the first and second derivatives of the friction dissipative potential with respect to the velocities.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::VectorXd friction_potential_grad =
                friction_constraints.compute_potential_gradient(
                    collision_mesh, velocity, epsv);

            Eigen::SparseMatrix<double> friction_potential_hess =
                friction_constraints.compute_potential_hessian(
                    collision_mesh, velocity, epsv);

    .. md-tab-item:: Python

        .. code-block:: python

            friction_potential_grad = friction_constraints.compute_potential_gradient(
                collision_mesh, velocity, epsv)

            friction_potential_hess = friction_constraints.compute_potential_hessian(
                collision_mesh, velocity, epsv)

Continuous Collision Detection
------------------------------

The last high-level component of the IPC Toolkit library is continuous collision detection (CCD). This is a method for determining if and at what time two objects will collide. This can be incorporated in a simulation nonlinear solver's line-search to determine the maximum step size allowable before a collision occurs.

There are two main functions for doing this: ``is_step_collision_free`` and ``compute_collision_free_stepsize``. The former returns a boolean value indicating if the step is collision-free, and the latter returns the maximum step size that is collision-free. Both functions take the same arguments, but ``compute_collision_free_stepsize`` is the more convenient function to use as it returns the maximum step size.

The following example determines the maximum step size allowable between the rest_positions and the squashed bunny.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::MatrixXd vertices_t0 = collision_mesh.rest_positions();
            Eigen::MatrixXd vertices_t1 = vertices_t0;
            vertices_t1.col(1) *= 0.01;

            double max_step_size = compute_collision_free_stepsize(
                    collision_mesh, vertices_t0, vertices_t1);

            Eigen::MatrixXd collision_free_vertices =
                (vertices_t1 - vertices_t0) * max_step_size + vertices_t0;
            assert(is_step_collision_free(mesh, vertices_t0, collision_free_vertices));

    .. md-tab-item:: Python

        .. code-block:: python

            vertices_t0 = collision_mesh.rest_positions()
            vertices_t1 = vertices_t0
            vertices_t1[:, 1] *= 0.01

            max_step_size = compute_collision_free_stepsize(
                    collision_mesh, vertices_t0, vertices_t1)

            collision_free_vertices =
                (vertices_t1 - vertices_t0) * max_step_size + vertices_t0
            assert(is_step_collision_free(mesh, vertices_t0, collision_free_vertices))