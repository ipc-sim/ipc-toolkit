Physical Simulation
===================

While the IPC Toolkit provides all the principle components of the IPC algorithm, it does not provide a complete simulation framework. Instead, it provides building blocks that can be used to integrate the IPC algorithm into a physical simulation. If all you want is a complete simulation framework using the IPC algorithm, then you should check out our other project `PolyFEM <https://polyfem.github.io/>`_ which uses the IPC Toolkit for its collision handling.

We provide several helper functions to make your job easier. The following examples show how to use these functions.

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
            Eigen::MatrixXd faces;
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
            Eigen::VectorXd vertices = collision_mesh.vertices(full_vertices);

            // Construct the set of collisions
            ipc::Collisions collisions;
            collisions.build(collision_mesh, vertices, dhat);

            // Construct a barrier potential
            ipc::BarrierPotential B(dhat);

            // Evaluate the potential
            double b = B(collisions, collision_mesh, vertices);

            // Convert full velocities to surface velocities
            Eigen::VectorXd velocities = collision_mesh.map_displacements(full_velocities);

            // Construct the set of friction collisions
            ipc::FrictionCollisions friction_collisions;
            friction_collisions.build(collision_mesh, vertices, collisions, B, barrier_stiffness, mu);

            // Construct a friction dissipative potential
            ipc::FrictionPotential D(epsv);

            double d = D(friction_collisions, collision_mesh, velocities);

    .. md-tab-item:: Python

        .. code-block:: python

            # Convert full vertices to surface vertices
            vertices = collision_mesh.vertices(full_vertices)

            # Construct the set of collisions
            collisions = ipctk.Collisions()
            collisions.build(collision_mesh, vertices, dhat)

            # Construct a barrier potential
            B = ipctk.BarrierPotential(dhat)

            # Evaluate the potential
            b = B(collisions, collision_mesh, vertices)

            # Convert full velocities to surface velocities
            velocities = collision_mesh.map_displacements(full_velocities)

            # Construct the set of friction collisions
            friction_collisions = ipctk.FrictionCollisions()
            friction_collisions.build(collision_mesh, vertices, collisions, B, barrier_stiffness, mu)

            # Construct a friction dissipative potential
            D = ipctk.FrictionPotential(epsv)

            d = D(friction_collisions, collision_mesh, velocities)

When computing the gradient and Hessian of the potentials, the derivatives will be with respect to the surface DOF. If you want the derivatives with respect to the full mesh DOF, then we need to apply the chain rule. Fortunately, the ``CollisionMesh`` class provides a function to do this (``CollisionMesh::to_full_dof``):

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            const BarrierPotential B(dhat);

            Eigen::VectorXd grad = B.gradient(collisions, collision_mesh, vertices);
            Eigen::VectorXd grad_full = collision_mesh.to_full_dof(grad);

            Eigen::SparseMatrix<double> hess = B.hessian(collisions, collision_mesh, vertices);
            Eigen::SparseMatrix<double> hess_full = collision_mesh.to_full_dof(hess);

    .. md-tab-item:: Python

        .. code-block:: python

            B = BarrierPotential(dhat)

            grad = B.gradient(collision, collision_mesh, vertices)
            grad_full = collision_mesh.to_full_dof(grad)

            hess = B.hessian(collision, collision_mesh, vertices)
            hess_full = collision_mesh.to_full_dof(hess)

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

            // Construct the collision mesh from the is_on_surface vector and full mesh data
            ipc::CollisionMesh collision_mesh(
                is_on_surface, full_rest_positions, edges, faces);

    .. md-tab-item:: Python

        .. code-block:: python

            # codim_vertices is an array of indices of the codimensional vertices
            codim_vertices = ...

            # is_on_surface is a list of booleans indicating which vertices are on the surface
            is_on_surface = ipctk.CollisionMesh.construct_is_on_surface(
                len(full_rest_positions), boundary_edges, codim_vertices)

            # Construct the collision mesh from the is_on_surface vector and full mesh data
            collision_mesh = ipctk.CollisionMesh(
                is_on_surface, full_rest_positions, edges, faces)

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
            igl::readOBJ("proxy.obj", rest_positions, faces);
            igl::edges(faces, edges);
            // Or build it from the volumetric mesh

            // Linear map from the finite element mesh to the collision proxy
            Eigen::SparseMatrix<double> displacement_map = ...; // build or load the displacement map

            ipc::CollisionMesh collision_mesh(
                proxy_rest_positions, proxy_edges, proxy_faces, displacement_map);

    .. md-tab-item:: Python

        .. code-block:: python

            # Finite element mesh
            fe_mesh = meshio.read("mesh.msh")
            fe_rest_positions = mesh.points
            tets = mesh.cells_dict["tetra"]

            # Collision proxy mesh
            # Load the proxy mesh from a file
            proxy_mesh = meshio.read("proxy.obj")
            proxy_rest_positions = proxy_mesh.points
            proxy_faces = proxy_mesh.cells_dict["triangle"]
            proxy_edges = igl.edges(proxy_faces)
            # or build it from the volumetric mesh ...

            # Linear map from the finite element mesh to the collision proxy
            displacement_map = ... # build or load the displacement map

            collision_mesh = CollisionMesh(
                proxy_rest_positions, proxy_edges, proxy_faces, displacement_map)

We can then map the displacements using ``collision_mesh.map_displacement(fe_displacements)`` or directly get the displaced proxy mesh vertices using ``collision_mesh.displace_vertices(fe_displacements)``. Similarly, we can map forces/potential gradients using ``collision_mesh.to_full_dof(collision_forces)`` or force Jacobians/potential Hessians using ``collision_mesh.to_full_dof(potential_hessian)``.

.. warning::
    The function ``CollisionMesh::vertices(full_positions)`` should not be used in this case because the rest positions used to construct the ``CollisionMesh`` are not the same as the finite element mesh's rest positions. Instead, use ``CollisionMesh::displace_vertices(fe_displacements)`` where ``fe_displacements`` is already the solution of the PDE or can be computed as ``fe_displacements = fe_positions - fe_rest_positions`` from deformed and rest positions.

Positive Semi-Definite Projection
---------------------------------

As described by :cite:t:`Li2020IPC`, the Hessian of the potentials can be indefinite. This is problematic when using the Hessian in a Newton step :cite:p:`Li2020IPC`.
To remedy this, we can project the Hessian onto the positive semidefinite (PSD) cone. To do this set the optional parameter ``project_hessian_to_psd`` in ``Potential::hessian`` to one of the following.

.. md-tab-set::

    .. md-tab-item:: C++

        - ``ProjectToPSD::CLAMP``: Clamp the negative eigenvalues of the Hessian to 0. This is the same as used by :cite:t:`Li2020IPC`.
        - ``ProjectToPSD::ABS``: Set the negative eigenvalues of the Hessian to their absolute value. This is the method proposed by :cite:t:`Chen2024Stabler`.

    .. md-tab-item:: Python

        - ``ProjectToPSD.CLAMP``: Clamp the negative eigenvalues of the Hessian to 0. This is the same as used by :cite:t:`Li2020IPC`.
        - ``ProjectToPSD.ABS``: Set the negative eigenvalues of the Hessian to their absolute value. This is the method proposed by :cite:t:`Chen2024Stabler`.