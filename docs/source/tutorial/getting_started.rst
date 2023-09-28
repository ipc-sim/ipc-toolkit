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
            faces = mesh.cells_dict["triangle"]
            edges = ipctk.edges(faces)

            collision_mesh = ipctk.CollisionMesh(rest_positions, edges, faces)

The matrix ``rest_positions`` contains the undeformed positions of the vertices. It should be :math:`|V| \times d` where :math:`|V|` is the number of vertices and :math:`d \in \{2, 3\}` is the dimension.
Each row of the ``edges`` and ``faces`` matrices contains the vertex IDs (row number in ``rest_positions``) of the edge or face's vertices.
The sizes of ``edges`` and ``faces`` are ``#E x 2`` and ``#F x 3`` respectively (``#E`` and ``#F`` are the number of edges and faces).

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

This returns a scalar value ``contact_potential`` which is the sum of the contact potentials for each active constraint.

Mathematically this is defined as

.. math::
   B(x) = \sum_{k \in C} b(d_k(x), \hat{d}),

where :math:`x` is our deformed vertex positions, :math:`C` is the active collision constraints, :math:`d_k` is the distance (squared) of the :math:`k`-th active constraint, and :math:`b` is IPC's C2-clamped log-barrier function.

.. note::
   This is **not** premultiplied by the barrier stiffness :math:`\kappa`.

Contact Potential Derivatives
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

To compute the adaptive barrier stiffness, we can use two functions: ``initial_barrier_stiffness`` and ``update_barrier_stiffness``. The function ``initial_barrier_stiffness``computes the initial value from the current energy and contact potential gradients. This function also provides a minimum and maximum value for the barrier stiffness. The function ``update_barrier_stiffness`` updates the barrier stiffness if the minimum distance has become too small.

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

.. _convergent-contact-formulation:

Convergent Formulation
^^^^^^^^^^^^^^^^^^^^^^

In addition to the original implementation of :cite:t:`Li2020IPC`, we also implement the convergent formulation of :cite:t:`Li2023Convergent`.

To enable the convergent formulation, we need to set ``use_convergent_formulation``:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            collision_constraints.set_use_convergent_formulation(true);
            collision_constraints.build(collision_mesh, vertices, dhat);

    .. md-tab-item:: Python

        .. code-block:: python

            collision_constraints.use_convergent_formulation = True
            collision_constraints.build(collision_mesh, vertices, dhat)

.. important::
    The variable ``use_convergent_formulation`` should be set before calling ``CollisionConstraints::build`` for it to take effect. By default, it is ``false``.

Technical Details
'''''''''''''''''

*We breifly summarize the convergent formulation here for convenience.*

In order to derive a convergent formulation, we first define a continous form of our barrier potential :math:`P`. For a surface :math:`\mathcal{S}` embeded in 3D space, we parameterize the surfaces by common (possibly discontinuous) coordinates :math:`u \in \tilde{M} \subset \mathbb{R}^2`, so that :math:`\mathcal{S}(u)` traverses the material points across all surfaces contiguously. The total contact potential is then

.. math::
   P(\mathcal{S})=\frac{1}{2} \int_{u \in \tilde{M}} \max _{v \in \tilde{M} \setminus{ }_r u} b(d(\mathcal{S}(u), \mathcal{S}(v)), \hat{d})~\mathrm{d} u,

where we define the operator :math:`\setminus_r: \mathcal{P}(\mathbb{R}^2) \times \mathbb{R} \times \mathbb{R}^2 \mapsto \mathcal{P}(\mathbb{R}^2)` to be

.. math::
    \tilde{M} \setminus_r u:=\left\{v \in \tilde{M} \mid\|u-v\|_2>r\right\}

with :math:`r \rightarrow 0`.

We then define our surface discretization with a triangulated boundary mesh geometry. As in the smooth case we can parameterize the domain across all polylines with :math:`u \in \tilde{M}` so that :math:`p(u): \tilde{M} \mapsto \mathbb{R}^3` traverses all material points, across all triangles :math:`t ∈ T` in the triangle mesh contiguously. The corresponding surface contact potential is then

.. math::
    \frac{1}{2} \int_{u \in \tilde{M}} \max_{t \in T \backslash p(u)} b(d(p(u), t), \hat{d})~\mathrm{d} u,

where :math:`T \setminus p` is the set of boundary faces that do not contain the point p.

Applying mesh vertices as nodes (and quadrature points), we numerically integrate the surface contact potential. For each nodal position :math:`x \in V` we then have a corresponding material space coordinate :math:`\bar{x} \in \bar{V}`. Piecewise linear integration of the surface barrier is then

.. math::
    \frac{1}{2} \sum_{\bar{x} \in \bar{V}} w_{\bar{x}} \max _{t \in T \backslash x(\bar{x})} b(d(x(\bar{x}), t), \hat{d}),

where :math:`w_{\bar{x}}` are the quadrature weights, each given by one third of the sum of the areas (in material space) of the boundary triangles incident to :math:`\bar{x}`.

We next need to smoothly approximate the max operator in the contact potentials. However, common approaches such as an :math:`L^p`-norm or LogSumExp would decrease sparsity in subsequent numerical solves by increasing stencil-size per contact evaluation. We instead leverage the locality of our barrier function to approximate the max operator by removing duplicate distance pairs. Our resulting approximators for a trianglated surfaces is

.. math::
    \begin{aligned}
    \Psi_s(x) & =\sum_{t \in T \backslash x} b(d(x, t), \hat{d})-\sum_{e \in E_{\text {int }} \backslash x} b(d(x, e), \hat{d})+\sum_{x_2 \in V_{i n t} \backslash x} b\left(d\left(x, x_2\right), \hat{d}\right) \\
    & \approx \max _{t \in T \backslash x} b(d(x, t), \hat{d}),
    \end{aligned}

where :math:`V_{\text{int}} \subseteq V` is the subset of internal surface nodes and :math:`E_{\text{int}} \subseteq E` is the subset of internal surface edges (i.e., edges incident to two triangles). For locally convex regions this estimator is tight while remaining smooth. In turn, for nonconvex regions it improves over direct summation.

The corresponding discrete barrier potential is then simply

.. math::
    P_s(V)= \frac{1}{2} \sum_{x \in V} w_x \Psi_s(x).

Please, see the `paper <https://arxiv.org/abs/2307.15908>`_ for more details (including the formulation for 2D curves and edge-edge contact) and evaluation.

The key difference between the original and the convergent formulations is that we (1) include area weights in the barrier potential and (2) include additional (negative) terms to cancel out the duplicate distance pairs. While this requires minor algorithmic changes, it produces considerably better results.

Physical Barrier
''''''''''''''''

We want our barrier potential to have the same units of as our elastic potential (e.g., :math:`\text{J}`). Together with the area weighting (discussed above) this means the barrier should has units of pressure times distance (e.g., :math:`\text{Pa} \cdot \text{m}`). That is,

.. math::
    \text{Pa} \cdot \text{m} \cdot \text{m}^2 = \frac{\text{N}}{\text{m}^2} \cdot \text{m} \cdot \text{m}^2 = \text{N} \cdot \text{m} = \text{J}.

To achieve this, (when using the convergent formulation) we modify the barrier function to have units of distance:

.. math::
    b(d, \hat{d})=\left\{\begin{array}{lr}
    -\hat{d}\left(\frac{d}{\hat{d}}-1\right)^2 \ln \left(\frac{d}{\hat{d}}\right), & 0<d<\hat{d} \\
    0 & d \geq \hat{d}
    \end{array}\right.

.. note::
    This is equivalent to the original barrier function of :cite:p:`Li2020IPC` times :math:`1/\hat{d}^3` when using squared distances. Therefore, to simplify the implementation we only implement the original barrier function and multiply all constraints by :math:`1/\hat{d}^3`.

The barrier stiffness (:math:`\kappa`) then has units of pressure (e.g., :math:`\text{Pa}`), the same as Young's modulus (:math:`E`) in elasticity.
This implies we can get good solver convergence even when using a fixed :math:`\kappa` by setting it relative to the material's Young's modulus (:math:`\kappa = 0.1 E` works well in many examples).
The intention is to treat the barrier as a thin elastic region around the mesh, and having consistent units makes it easier to pick the stiffness for this "material".

.. _modeling-thickness:

Modeling Thickness
^^^^^^^^^^^^^^^^^^

We implement the thickness model of :cite:t:`Li2021CIPC` to apply an offset (refered to as :math:`\xi` in :cite:p:`Li2021CIPC` or :math:`d_\min` here) to the collision constraints. This is useful for modeling the thickness of a shell or cloth.

To add a constraint offset, we need to set the ``dmin`` variable. For example, we can set the constraint offset :math:`d_\min=10^{-3}` and :math:`\hat{d}=10^{-4}`:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            const double dhat = 1e-4;
            const double dmin = 1e-3;

            ipc::CollisionConstraints collision_constraints;
            collision_constraints.build(collision_mesh, vertices, dhat, dmin);

    .. md-tab-item:: Python

        .. code-block:: python

            dhat = 1e-4
            dmin = 1e-3

            collision_constraints = ipctk.CollisionConstraints()
            collision_constraints.build(collision_mesh, vertices, dhat, dmin)

This will then set the ``dmin`` field in all of the ``CollisionConstraint`` objects stored in the ``collision_constraints``.

.. note::
    Currently, only a single thickness value is supported for the entire mesh.

It is also important to set use the same :math:`d_\min` when performing CCD (see :ref:`Minimum Separation CCD <minimum-separation-ccd>`).

Friction
--------

Computing the friction dissipative potential is similar to the contact potential, but because it is a lagged model, we need to build it from a fixed set of constraints.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            ipc::FrictionConstraints friction_constraints;
            friction_constraints.build(
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
      v = \frac{x - x^t}{h},

   where :math:`x` is the current position, :math:`x^t` is the previous position, and :math:`h` is the time step size.

This returns a scalar value ``friction_potential`` which is the sum of the individual friction potentials.

Mathematically this is defined as

.. math::
   D(x) = \sum_{k \in C} \mu\lambda_k^nf_0\left(\|T_k^Tv\|, \epsilon_v\right),

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

.. _convergent-friction-formulation:

Convergent Formulation
^^^^^^^^^^^^^^^^^^^^^^

Just as with the :ref:`collision constraints <convergent-contact-formulation>`, we implement both the original friction formulation of :cite:t:`Li2020IPC` and the convergent formulation of :cite:t:`Li2023Convergent`.

The choice of formulation is dependent on how the fixed set of ``contact_constraints`` given to ``FrictionConstraints::build`` was built. If the ``contact_constraints`` were built using the convergent formulation, then the friction constraints will also use the convergent formulation. Otherwise, the original formulation will be used.

Continuous Collision Detection
------------------------------

The last high-level component of the IPC Toolkit library is continuous collision detection (CCD). This is a method for determining if and at what time two objects will collide. This can be incorporated in a simulation nonlinear solver's line-search to determine the maximum step size allowable before a collision occurs.

There are two main functions for doing this: ``is_step_collision_free`` and ``compute_collision_free_stepsize``. The former returns a boolean value indicating if the step is collision-free, and the latter returns the maximum step size that is collision-free. Both functions take the same arguments, but ``compute_collision_free_stepsize`` is the more convenient function to use as it returns the maximum step size.

The following example determines the maximum step size allowable between the rest_positions and the squashed bunny.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::MatrixXd vertices_t0 = collision_mesh.rest_positions(); // vertices at t=0
            Eigen::MatrixXd vertices_t1 = vertices_t0;                     // vertices at t=1
            vertices_t1.col(1) *= 0.01; // squash the mesh in the y-direction

            double max_step_size = compute_collision_free_stepsize(
                    collision_mesh, vertices_t0, vertices_t1);

            Eigen::MatrixXd collision_free_vertices =
                (vertices_t1 - vertices_t0) * max_step_size + vertices_t0;
            assert(is_step_collision_free(mesh, vertices_t0, collision_free_vertices));

    .. md-tab-item:: Python

        .. code-block:: python

            vertices_t0 = collision_mesh.rest_positions() # vertices at t=0
            vertices_t1 = vertices_t0.copy()              # vertices at t=1
            vertices_t1[:, 1] *= 0.01 # squash the mesh in the y-direction

            max_step_size = compute_collision_free_stepsize(
                    collision_mesh, vertices_t0, vertices_t1)

            collision_free_vertices =
                (vertices_t1 - vertices_t0) * max_step_size + vertices_t0
            assert(is_step_collision_free(mesh, vertices_t0, collision_free_vertices))

CCD is comprised of two parts (phases): broad-phase and narrow-phase.

Broad-Phase
^^^^^^^^^^^

The broad-phase takes all possible pairings (quadratic in size) and eliminates (culls) pairs whose bounding boxes do not overlap. This is done using a spatial data structure (e.g., a hash grid or spatial hash).

The ``Candidates`` class represents the culled set of candidate pairs and is built by using a broad-phase method. The following example shows how to use the broad-phase to determine the candidate pairs between the rest_positions and the squashed bunny.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <ipc/candidates/candidates.hpp>

            ipc::Candidates candidates;
            candidates.build(
                mesh, vertices_t0, vertices_t1,
                /*inflation_radius=*/0.0,
                /*broad_phase_method=*/ipc::BroadPhaseMethod::HASH_GRID);

    .. md-tab-item:: Python

        .. code-block:: python

            candidates = ipctk.Candidates()
            candidates.build(
                mesh, vertices_t0, vertices_t1,
                broad_phase_method=ipctk.BroadPhaseMethod.HASH_GRID)

Possible values for ``broad_phase_method`` are: ``BRUTE_FORCE`` (parallel brute force culling), ``HASH_GRID`` (default), ``SPATIAL_HASH`` (implementation from the original IPC codebase), ``SWEEP_AND_TINIEST_QUEUE`` (method of :cite:t:`Belgrod2023Time`), or ``SWEEP_AND_TINIEST_QUEUE_GPU`` (requires CUDA).

Narrow-Phase
^^^^^^^^^^^^

The narrow-phase computes the time of impact between two primitives (e.g., a point and a triangle or two edges in 3D). To do this we utilize the Tight Inclusion CCD method of :cite:t:`Wang2021TightInclusion` for the narrow-phase as it is provably conservative (i.e., never misses collisions), accurate (i.e., rarely reports false positives), and efficient.

The following example shows how to use the narrow-phase to determine if a point is colliding with a triangle (static in this case).

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            #include <ipc/ccd/ccd.hpp>

            // ...

            Eigen::Vector3d p_t0(0.0, -1.0, 0.0); // point at t=0
            Eigen::Vector3d p_t1(0.0,  1.0, 0.0); // point at t=1

            Eigen::Vector3d t0_t0(-1.0, 0.0,  1.0); // triangle vertex 0 at t=0
            Eigen::Vector3d t1_t0( 1.0, 0.0,  1.0); // triangle vertex 1 at t=0
            Eigen::Vector3d t2_t0( 0.0, 0.0, -1.0); // triangle vertex 2 at t=0

            // static triangle
            Eigen::Vector3d t0_t1 = t0_t0; // triangle vertex 0 at t=1
            Eigen::Vector3d t1_t1 = t1_t0; // triangle vertex 1 at t=1
            Eigen::Vector3d t2_t1 = t2_t0; // triangle vertex 2 at t=1

            double toi; // output time of impact
            bool is_colliding = ipc::point_triangle_ccd(
                p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, toi);
            assert(is_colliding);
            assert(abs(toi - 0.5) < 1e-8);

    .. md-tab-item:: Python

        .. code-block:: python

            import numpy as np
            import ipctk

            p_t0 = np.array([0.0, -1.0, 0.0]) # point at t=0
            p_t1 = np.array([0.0,  1.0, 0.0]) # point at t=1

            t0_t0 = np.array([-1.0, 0.0,  1.0]) # triangle vertex 0 at t=0
            t1_t0 = np.array([ 1.0, 0.0,  1.0]) # triangle vertex 1 at t=0
            t2_t0 = np.array([ 0.0, 0.0, -1.0]) # triangle vertex 2 at t=0

            # static triangle
            t0_t1 = t0_t0 # triangle vertex 0 at t=1
            t1_t1 = t1_t0 # triangle vertex 1 at t=1
            t2_t1 = t2_t0 # triangle vertex 2 at t=1

            # returns a boolean indicating if the point is colliding with the triangle
            # and the time of impact (TOI)
            is_colliding, toi = ipctk.point_triangle_ccd(
                p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1)
            assert(is_colliding)
            assert(abs(toi - 0.5) < 1e-8)

Alternatively, the ``FaceVertexCandidate`` class contains a ``ccd`` function that can be used to determine if the face-vertex pairing are colliding:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            ipc::FaceVertexCandidate candidate = ...; // face-vertex candidate

            double toi; // output time of impact
            bool is_colliding = candidate.ccd(
                vertices_t0, vertices_t1, collision_mesh.edges(), collision_mesh.faces(), toi);

    .. md-tab-item:: Python

        .. code-block:: python

            candidate = ... # face-vertex candidate

            # returns a boolean indicating if the point is colliding with the triangle
            # and the time of impact (TOI)
            is_colliding, toi = candidate.ccd(
                vertices_t0, vertices_t1, collision_mesh.edges, collision_mesh.faces)

The same can be done for point-edge collisions using the ``point_edge_ccd`` function or ``EdgeVertexCandidate`` class and for edge-edge collisions using the ``edge_edge_ccd`` function or ``EdgeEdgeCandidate`` class.

.. _minimum-separation-ccd:

Minimum Separation
^^^^^^^^^^^^^^^^^^

We can also perform CCD with a minimum separation distance. This is useful when modeling thickness (see, e.g., :ref:`Modeling Thickness <modeling-thickness>`).

To do this, we need to set the ``min_distance`` parameter when calling ``is_step_collision_free`` and ``compute_collision_free_stepsize``. For example, we can set the minimum separation distance to :math:`d_\min=10^{-4}`:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            double max_step_size = compute_collision_free_stepsize(
                    collision_mesh, vertices_t0, vertices_t1,
                    /*broad_phase_method=*/ipc::DEFAULT_BROAD_PHASE_METHOD,
                    /*min_distance=*/1e-4);

            Eigen::MatrixXd collision_free_vertices =
                (vertices_t1 - vertices_t0) * max_step_size + vertices_t0;
            assert(is_step_collision_free(
                mesh, vertices_t0, collision_free_vertices,
                /*broad_phase_method=*/ipc::DEFAULT_BROAD_PHASE_METHOD,
                /*min_distance=*/1e-4
            ));

    .. md-tab-item:: Python

        .. code-block:: python

            max_step_size = compute_collision_free_stepsize(
                    collision_mesh, vertices_t0, vertices_t1, min_distance=1e-4)

            collision_free_vertices =
                (vertices_t1 - vertices_t0) * max_step_size + vertices_t0
            assert(is_step_collision_free(
                mesh, vertices_t0, collision_free_vertices, min_distance=1e-4))