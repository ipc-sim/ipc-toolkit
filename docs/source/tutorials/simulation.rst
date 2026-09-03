Physical Simulation
===================

While the IPC Toolkit provides all the principle components of the IPC algorithm, it does not provide a complete simulation framework. Instead, it provides building blocks that can be used to integrate the IPC algorithm into a physical simulation. If all you want is a complete simulation framework using the IPC algorithm, then you should check out our other project `PolyFEM <https://polyfem.github.io/>`_ which uses the IPC Toolkit for its collision handling.

We provide several helper functions to make your job easier. The following examples show how to use these functions.

Optimization-Based Time Integration
-----------------------------------

IPC defines barrier potential :math:`B(\mathbf{x})` and a friction potential :math:`D(\mathbf{v})`. To add these into a optimization-based time integration, we need to scale the potentials by the time-integrators acceleration scaling. For implicit Euler, this is h^2, where h is the timestep.

With elasticity :math:`\Psi(\mathbf{x})`, the total optimization problem is:

.. math::
   \mathbf{x}^{t+1} = \underset{\mathbf{x}}{\arg\min} ~ \tfrac{1}{2} (\mathbf{x} - \hat{\mathbf{x}})^\top\mathbf{M}(\mathbf{x}-\hat{\mathbf{x}})+h^2\Psi(\mathbf{x}) + h^2 \kappa B(\mathbf{x}) + h^2 D(\mathbf{v}(\mathbf{x}))

where :math:`\hat{\mathbf{x}} = \mathbf{x}^t + h\mathbf{v}^t + h^2\mathbf{g}` is the time integration scheme-specific “predicted positions.”

.. note::
    In :cite:t:`Li2020IPC`, all the constants are wrapped up into $\kappa$, which is adaptively modified. In follow-up works, we treat the barrier as a physical energy, and so it should have the same multiplier as the elastic energy ($h^2$ for implicit Euler).

Volumetric Meshes
-----------------

The IPC Toolkit only handles surface meshes (through the ``CollisionMesh``). However, the finite element method often relies on volumetric discretization. In this case, the computed gradients and Hessians need to be mapped back to the full volumetric mesh. The ``CollisionMesh`` class provides this functionality.

From the full (volumetric) mesh vertices and surface edges/faces which index into the full mesh vertices, you can build a ``CollisionMesh`` using the function ``CollisionMesh::build_from_full_mesh``. This will internally build and store a selection matrix that goes from the full to surface vertices as well as map the edge/faces entries accordingly.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::MatrixXd full_rest_positions;
            Eigen::MatrixXi tets;
            // TODO: Show how to load a volumetric mesh from a file (e.g., using MshIO)

            // Faces of the surface mesh with indices into full_rest_positions
            Eigen::MatrixXi faces;
            igl::boundary_facets(tets, faces);

            // Edges of the surface mesh with indices into full_rest_positions
            Eigen::MatrixXi edges;
            igl::edges(faces, edges);

            ipc::CollisionMesh collision_mesh =
                ipc::CollisionMesh::build_from_full_mesh(full_rest_positions, edges, faces);

    .. md-tab-item:: Python

        .. code-block:: python

            mesh = meshio.read("bunny.msh")
            full_rest_positions = mesh.points
            tets = mesh.cells_dict["tetra"]

            faces = igl.boundary_facets(tets)  # pip install libigl
            edges = ipctk.edges(faces)         # same as igl.edges

            collision_mesh = ipctk.CollisionMesh.build_from_full_mesh(
                full_rest_positions, edges, faces)

This ``CollisionMesh`` can then be used just as any other ``CollisionMesh``. However, when passing the collision mesh to toolkit functions, the vertices have to be the surface vertices. The ``CollisionMesh`` class provides a function to map the full vertices and velocities to the surface vertices and velocities. These are ``CollisionMesh::vertices(full_vertices)`` and ``CollisionMesh::map_displacements(full_displacements)``, respectively.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // Convert full vertices to surface vertices
            Eigen::MatrixXd vertices = collision_mesh.vertices(full_vertices);

            // Construct the set of collisions
            ipc::NormalCollisions collisions;
            collisions.build(collision_mesh, vertices, dhat);

            // Construct a barrier potential
            ipc::BarrierPotential B(dhat, stiffness);

            // Evaluate the potential
            double b = B(collisions, collision_mesh, vertices);

            // Convert full velocities to surface velocities
            Eigen::MatrixXd velocities = collision_mesh.map_displacements(full_velocities);

            // Construct the set of friction collisions
            ipc::TangentialCollisions tangential_collisions;
            tangential_collisions.build(collision_mesh, vertices, collisions, B, mu);

            // Construct a friction dissipative potential
            ipc::FrictionPotential D(eps_v);

            double d = D(tangential_collisions, collision_mesh, velocities);

    .. md-tab-item:: Python

        .. code-block:: python

            # Convert full vertices to surface vertices
            vertices = collision_mesh.vertices(full_vertices)

            # Construct the set of collisions
            collisions = ipctk.NormalCollisions()
            collisions.build(collision_mesh, vertices, dhat)

            # Construct a barrier potential
            B = ipctk.BarrierPotential(dhat, stiffness)

            # Evaluate the potential
            b = B(collisions, collision_mesh, vertices)

            # Convert full velocities to surface velocities
            velocities = collision_mesh.map_displacements(full_velocities)

            # Construct the set of friction collisions
            tangential_collisions = ipctk.TangentialCollisions()
            tangential_collisions.build(collision_mesh, vertices, collisions, B, mu)

            # Construct a friction dissipative potential
            D = ipctk.FrictionPotential(eps_v)

            d = D(tangential_collisions, collision_mesh, velocities)

When computing the gradient and Hessian of the potentials, the derivatives will be with respect to the surface DOF. If you want the derivatives with respect to the full mesh DOF, then we need to apply the chain rule. Fortunately, the ``CollisionMesh`` class provides a function to do this (``CollisionMesh::to_full_dof``):

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            const BarrierPotential B(dhat, stiffness);

            Eigen::VectorXd grad = B.gradient(collisions, collision_mesh, vertices);
            Eigen::VectorXd grad_full = collision_mesh.to_full_dof(grad);

            Eigen::SparseMatrix<double> hess = B.hessian(collisions, collision_mesh, vertices);
            Eigen::SparseMatrix<double> hess_full = collision_mesh.to_full_dof(hess);

    .. md-tab-item:: Python

        .. code-block:: python

            B = ipctk.BarrierPotential(dhat, stiffness)

            grad = B.gradient(collisions, collision_mesh, vertices)
            grad_full = collision_mesh.to_full_dof(grad)

            hess = B.hessian(collisions, collision_mesh, vertices)
            hess_full = collision_mesh.to_full_dof(hess)

If only the full DOF derivatives are required, the optional ``in_full_dof`` parameter requests them directly rather than mapping after the fact. The stencil indices are remapped during assembly, so the gradient is scattered into full DOF as it is accumulated and the Hessian avoids the two sparse matrix products performed by ``to_full_dof``:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::VectorXd grad_full = B.gradient(
                collisions, collision_mesh, vertices, /*in_full_dof=*/true);

            Eigen::SparseMatrix<double> hess_full = B.hessian(
                collisions, collision_mesh, vertices,
                ipc::PSDProjectionMethod::NONE, /*in_full_dof=*/true);

    .. md-tab-item:: Python

        .. code-block:: python

            grad_full = B.gradient(
                collisions, collision_mesh, vertices, in_full_dof=True)

            hess_full = B.hessian(
                collisions, collision_mesh, vertices, in_full_dof=True)

The results agree with ``collision_mesh.to_full_dof(...)`` up to the order in which the local contributions are summed.

.. note::
    Remapping the indices is only valid when the map from collision to full DOF is a pure selection, which holds for any collision mesh constructed without a displacement map. Use ``collision_mesh.is_selection_dof_map()`` to query this. When a displacement map is present, ``in_full_dof`` still returns the correct result, but it does so by applying ``to_full_dof`` internally and therefore offers no advantage.

Codimensional Vertices
^^^^^^^^^^^^^^^^^^^^^^

In some cases, the collision mesh vertices are not the same as the surface vertices of the volumetric mesh vertices. One such case is when simulating codimensional vertices in conjunction with shell or volumetric meshes. In this case, simply calling ``build_from_full_mesh`` will not work as it will ignore the vertices that are not connected to any boundary edge. Instead, you can build a vector of booleans that indicate which vertices are on the surface and pass it to the ``CollisionMesh`` constructor.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // codim_vertices is a vector of indices of the codimensional vertices
            Eigen::VectorXi codim_vertices = ...;

            // is_on_surface is a vector of booleans indicating which vertices are on the surface
            std::vector<bool> is_on_surface = ipc::CollisionMesh::construct_is_on_surface(
                full_rest_positions.rows(), boundary_edges, codim_vertices);

            // is_orient_vertex marks the vertices with a well-defined orientation
            // (i.e., those usable for signed distances). Codimensional vertices
            // have no orientation, so false is a safe default for all vertices.
            std::vector<bool> is_orient_vertex(full_rest_positions.rows(), false);

            // Construct the collision mesh from the masks and full mesh data
            ipc::CollisionMesh collision_mesh(
                is_on_surface, is_orient_vertex, full_rest_positions, edges, faces);

    .. md-tab-item:: Python

        .. code-block:: python

            # codim_vertices is an array of indices of the codimensional vertices
            codim_vertices = ...

            # is_on_surface is a list of booleans indicating which vertices are on the surface
            is_on_surface = ipctk.CollisionMesh.construct_is_on_surface(
                len(full_rest_positions), boundary_edges, codim_vertices)

            # is_orient_vertex marks the vertices with a well-defined orientation
            # (i.e., those usable for signed distances). Codimensional vertices
            # have no orientation, so False is a safe default for all vertices.
            is_orient_vertex = [False] * len(full_rest_positions)

            # Construct the collision mesh from the masks and full mesh data
            collision_mesh = ipctk.CollisionMesh(
                is_on_surface, is_orient_vertex, full_rest_positions, edges, faces)

Nonlinear Bases and Curved Meshes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

While IPC cannot directly handle nonlinear finite element bases and/or curved meshes, :cite:t:`Ferguson2023HighOrderIPC` show that displacements and forces can be transferred between a finite element mesh and a collision proxy through the use of a linear map. Given this linear map as a matrix, we can use the ``CollisionMesh`` class to map between the full and surface DOF.

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // Finite element mesh
            Eigen::MatrixXd fe_rest_positions;
            Eigen::MatrixXi tets;
            // TODO: Show how to load a volumetric mesh from a file (e.g., using MshIO)

            // Collision proxy mesh
            Eigen::MatrixXd proxy_rest_positions;
            Eigen::MatrixXi proxy_edges, proxy_faces;
            // Load the proxy mesh from a file
            igl::read_triangle_mesh("proxy.ply", proxy_rest_positions, proxy_faces);
            igl::edges(proxy_faces, proxy_edges);
            // Or build it from the volumetric mesh

            // Linear map from the finite element mesh to the collision proxy
            Eigen::SparseMatrix<double> displacement_map = ...; // build or load the displacement map

            ipc::CollisionMesh collision_mesh(
                proxy_rest_positions, proxy_edges, proxy_faces, displacement_map);

    .. md-tab-item:: Python

        .. code-block:: python

            # Finite element mesh
            fe_mesh = meshio.read("mesh.msh")
            fe_rest_positions = fe_mesh.points
            tets = fe_mesh.cells_dict["tetra"]

            # Collision proxy mesh
            # Load the proxy mesh from a file
            proxy_mesh = meshio.read("proxy.ply")
            proxy_rest_positions = proxy_mesh.points
            proxy_faces = proxy_mesh.cells_dict["triangle"]
            proxy_edges = igl.edges(proxy_faces)
            # or build it from the volumetric mesh ...

            # Linear map from the finite element mesh to the collision proxy
            displacement_map = ... # build or load the displacement map

            collision_mesh = ipctk.CollisionMesh(
                proxy_rest_positions, proxy_edges, proxy_faces, displacement_map)

We can then map the displacements using ``collision_mesh.map_displacements(fe_displacements)`` or directly get the displaced proxy mesh vertices using ``collision_mesh.displace_vertices(fe_displacements)``. Similarly, we can map forces/potential gradients using ``collision_mesh.to_full_dof(collision_forces)`` or force Jacobians/potential Hessians using ``collision_mesh.to_full_dof(potential_hessian)``.

.. warning::
    The function ``CollisionMesh::vertices(full_positions)`` should not be used in this case because the rest positions used to construct the ``CollisionMesh`` are not the same as the finite element mesh's rest positions. Instead, use ``CollisionMesh::displace_vertices(fe_displacements)`` where ``fe_displacements`` is already the solution of the PDE or can be computed as ``fe_displacements = fe_positions - fe_rest_positions`` from deformed and rest positions.

Positive Semi-Definite Projection
---------------------------------

As described by :cite:t:`Li2020IPC`, the Hessian of the potentials can be indefinite. This is problematic when using the Hessian in a Newton step :cite:p:`Li2020IPC`.
To remedy this, we can project the Hessian onto the positive semidefinite (PSD) cone. To do this set the optional parameter ``project_hessian_to_psd`` in ``Potential::hessian`` to one of the following.

.. md-tab-set::

    .. md-tab-item:: C++

        - ``PSDProjectionMethod::NONE``: Do not project the Hessian. This is the default.
        - ``PSDProjectionMethod::CLAMP``: Clamp the negative eigenvalues of the Hessian to 0. This is the same as used by :cite:t:`Li2020IPC`.
        - ``PSDProjectionMethod::ABS``: Set the negative eigenvalues of the Hessian to their absolute value. This is the method proposed by :cite:t:`Chen2024Stabler`.

        .. code-block:: c++

            Eigen::SparseMatrix<double> hess = B.hessian(
                collisions, collision_mesh, vertices,
                ipc::PSDProjectionMethod::CLAMP);

    .. md-tab-item:: Python

        - ``PSDProjectionMethod.NONE``: Do not project the Hessian. This is the default.
        - ``PSDProjectionMethod.CLAMP``: Clamp the negative eigenvalues of the Hessian to 0. This is the same as used by :cite:t:`Li2020IPC`.
        - ``PSDProjectionMethod.ABS``: Set the negative eigenvalues of the Hessian to their absolute value. This is the method proposed by :cite:t:`Chen2024Stabler`.

        .. code-block:: python

            hess = B.hessian(
                collisions, collision_mesh, vertices,
                project_hessian_to_psd=ipctk.PSDProjectionMethod.CLAMP)

Reusing the Hessian Assembler
-----------------------------

Each call to ``Potential::hessian`` constructs a sparse matrix from scratch. Except on small scenes, evaluating the local Hessians accounts for a minority of the cost; the bulk is spent determining where each local contribution belongs in the global matrix. A Newton solve repeats that work every iteration, even though the contact set typically changes little between iterations.

``Potential::assemble_hessian`` accepts the assembler as a parameter, allowing a single instance to persist across the solve and retain its sparsity pattern:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            // A single assembler for the entire solve, rather than one per iteration.
            ipc::MeshFEMHessianAssembler assembler;

            for (int i = 0; i < max_iterations; i++) {
                // ... update vertices and rebuild the collision set ...

                B.assemble_hessian(
                    collisions, collision_mesh, vertices, assembler,
                    ipc::PSDProjectionMethod::CLAMP, /*in_full_dof=*/true);

                // Valid until the next assembly; copy it to retain it longer.
                const Eigen::SparseMatrix<double>& hess = assembler.get_matrix();

                // ... solve for the Newton direction, line search, etc. ...
            }

    .. md-tab-item:: Python

        .. code-block:: python

            # A single assembler for the entire solve, rather than one per iteration.
            assembler = ipctk.MeshFEMHessianAssembler()

            for i in range(max_iterations):
                # ... update vertices and rebuild the collision set ...

                B.assemble_hessian(
                    collisions, collision_mesh, vertices, assembler,
                    ipctk.PSDProjectionMethod.CLAMP, in_full_dof=True)

                hess = assembler.get_matrix()

                # ... solve for the Newton direction, line search, etc. ...

The first call traverses the collision stencils and builds a block sparsity pattern, allocating one :math:`d \times d` block per interacting vertex pair instead of one entry per scalar. Subsequent calls compare the new stencils against the cached pattern. A newly active contact introduces a block the pattern does not contain and therefore forces a rebuild. A separating contact merely leaves behind a block that assembles to zero, which consumes some memory but does not alter the matrix, so the pattern is retained provided no more than ``stale_block_tolerance`` blocks have become stale.

The assembled matrix is symmetric and stored with only its upper triangle; ``get_matrix()`` mirrors it into a full Eigen matrix, reusing the cached structure whenever the pattern is unchanged. Solvers that consume block CSC directly can instead call ``block_matrix()`` to obtain MeshFEM's representation and avoid the conversion. Because our header only forward declares that type, such callers must include ``<MeshFEMSparse/BlockCSCHessian.hh>`` themselves.

When reassembling without modifying the collision set, for example under a different stiffness or PSD projection, the comparison itself can also be skipped:

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            assembler.set_assume_unchanged_stencils(true);

    .. md-tab-item:: Python

        .. code-block:: python

            assembler.assume_unchanged_stencils = True

.. warning::
    ``assume_unchanged_stencils`` is an unchecked assertion. If the stencils did change while their count remained equal, assembly reads past the end of the pattern. Debug builds verify the assumption; release builds do not.

.. note::
    ``MeshFEMHessianAssembler`` requires the toolkit to be built with ``IPC_TOOLKIT_WITH_MESHFEM_SPARSE`` (enabled by default, and the backend ``Potential::hessian`` uses internally). ``TripletHessianAssembler`` is always available and reproduces the historical behavior of ``hessian``, accumulating thread-local triplets and merging them with ``setFromTriplets``. It maintains no sparsity pattern, so reusing an instance of it confers no benefit.
