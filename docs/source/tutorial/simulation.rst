Physical Simulation
===================

While the IPC Toolkit provides all the principle components of the IPC algorithm, it does not provide a complete simulation framework. Instead, it provides building blocks that can be used to integrate the IPC algorithm into a physical simulation. If all you want is a complete simulation framework using the IPC algorithm, then you should check out our other project `PolyFEM <https://polyfem.github.io/>`_ which uses the IPC Toolkit for its collision handling.

We provide several helper functions to make your job easier. The following examples show how to use these functions.

Volumetric Meshes
-----------------

The IPC Toolkit only handles surface meshes (through the ``CollisionMesh``). However, the finite element method often relies on volumetric discretizations. In this case the computed gradients and Hessians need to be mapped back to the full volumetric mesh. The ``CollisionMesh`` class provides this functionality.

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

            std::vector<bool> is_on_surface = ipc::CollisionMesh::construct_is_on_surface(node_positions.rows(), boundary_edges);

            ipc::CollisionMesh collision_mesh =
                CollisionMesh::build_from_full_mesh(full_rest_positions, edges, faces);

    .. md-tab-item:: Python

        .. code-block:: python

            mesh = meshio.read("bunny.msh")
            full_rest_positions = mesh.points
            tets = mesh.cells[0].data

            faces = igl.boundary_facets(tets)  # pip install libigl
            edges = ipctk.edges(faces)         # same as igl.edges

            collision_mesh = ipctk.CollisionMesh.build_from_full_mesh(
                rest_positions, edges, faces)

This ``CollisionMesh`` can then be used just as any other ``CollisionMesh``. However, when computing the gradient and Hessian of the potentials, the derivatives will be with respect to the surface DOF. If you want the derivatives with respect to the full mesh DOF, then we need to apply the chain rule. Fortunately, the ``CollisionMesh`` class provides a function to do this (``CollisionMesh::to_full_dof``):

.. md-tab-set::

    .. md-tab-item:: C++

        .. code-block:: c++

            Eigen::VectorXd grad = collision_constraints.compute_potential_gradient(
                collision_mesh, vertices, dhat);
            Eigen::VectorXd grad_full = collision_mesh.to_full_dof(grad);

            Eigen::SparseMatrix<double> hess = collision_constraints.compute_potential_hessian(
                collision_mesh, vertices, dhat);
            Eigen::SparseMatrix<double> hess_full = collision_mesh.to_full_dof(hess);

    .. md-tab-item:: Python

        .. code-block:: python

            grad = collision_constraints.compute_potential_gradient(
                collision_mesh, vertices, dhat);
            grad_full = collision_mesh.to_full_dof(grad);

            hess = collision_constraints.compute_potential_hessian(
                collision_mesh, vertices, dhat);
            hess_full = collision_mesh.to_full_dof(hess);

Positive Semi-Definite Projection
---------------------------------

As described in [IPC]_ the Hessian of the potentials can be indefinite. This is problematic when using the Hessian in a Newton step [IPC]_. To remedy this, we can project the Hessian onto the positive semidefinite (PSD) cone. To do this set the optional parameter ``project_hessian_to_psd`` of ``compute_potential_hessian`` to true.

------------

.. rubric:: References

.. [IPC] Minchen Li, Zachary Ferguson, Teseo Schneider, Timothy Langlois, Denis Zorin, Daniele Panozzo, Chenfanfu Jiang, Danny M. Kaufman. 2020. Incremental Potential Contact: Intersection- and Inversion-free Large Deformation Dynamics. *ACM Transactions on Graphics (SIGGRAPH).*