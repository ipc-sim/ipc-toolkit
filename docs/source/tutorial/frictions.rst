Multiple Friction Materials
===========================

This tutorial demonstrates how to set up and work with multiple friction types in the IPC Toolkit. It covers setting up a defining tangential collisions with different friction properties, specifying blend types for friction coefficients, and calculating friction forces.


Create a Collision Mesh with Materials
-----------------------------------------

First, we need to create a collision mesh. The ``CollisionMesh`` data structure represents the surface geometry used for collision processing.
We will start by creating a collision mesh from a ``bunny.ply`` mesh file (you can find the mesh `here <https://github.com/ipc-sim/ipc-toolkit-tests-data/blob/main/bunny.ply>`_):

Adding a material ID to the ``CollisionMesh`` allows us to define friction properties for different materials.
The material ID is an integer value used to identify the material of each vertex in the mesh. It is recommended to use consecutive integers starting from 1 and
-1 for vertices that do not belong to any material.


.. code-block:: cpp

    #include <ipc/ipc.hpp>
    #include <Eigen/Core>
    #include <igl/read_triangle_mesh.h>
    #include <igl/edges.h>

    Eigen::MatrixXd rest_positions;
    Eigen::MatrixXi edges, faces;
    igl::read_triangle_mesh("bunny.ply", rest_positions, faces);
    igl::edges(faces, edges);

    // Assign material IDs to vertices
    std::vector<int> mat_id(rest_positions.rows(), 1);

    ipc::CollisionMesh collision_mesh(rest_positions, edges, faces);

In Python:

.. code-block:: python

    import ipctk
    import meshio

    mesh = meshio.read("bunny.ply")
    rest_positions = mesh.points
    faces = mesh.cells_dict["triangle"]
    edges = ipctk.edges(faces)

    collision_mesh = ipctk.CollisionMesh(rest_positions, edges, faces, mat_id=1)


Define Tangential Collisions with Multiple Friction Types and Blending
--------------------------------------------------------------------------

Define and build the tangential collisions for handling friction based on material properties. Specify `BlendType` to set the method for combining friction coefficients of different materials.

Supported blend types:
- `AVG`: Average
- `MIN`: Minimum
- `MAX`: Maximum
- `PRODUCT`: Product
- `HARMONIC_MEAN`: Harmonic mean
- `GEOMETRIC_MEAN`: Geometric mean

.. code-block:: cpp

    const double dhat = 1e-3;
    const double barrier_stiffness = 1.0;
    const double global_mu = 0.6;

    ipc::NormalCollisions collisions;
    collisions.build(collision_mesh, vertices, dhat);

    // Choose a blend type
    ipc::TangentialCollisions::BlendType blend_type = ipc::TangentialCollisions::BlendType::AVG;

    // blend function
    double default_blend_mu(double mu0, double mu1, ipc::TangentialCollisions::BlendType::AVG type)
    {
        switch (type) {
            case BlendType::MIN:
                return std::min(mu0, mu1);
            case BlendType::MAX:
                return std::max(mu0, mu1);
            case BlendType::PRODUCT:
                return mu0 * mu1;
            case BlendType::HARMONIC_MEAN:
                return 2 * (mu0 * mu1) / (mu0 + mu1);
            case BlendType::GEOMETRIC_MEAN:
                return std::sqrt(mu0 * mu1);
            case BlendType::AVG:
            default:
                return (mu0 + mu1) / 2;
        }
    }

    ipc::TangentialCollisions tangential_collisions;
    tangential_collisions.build(
        collision_mesh, vertices, collisions, ipc::BarrierPotential(dhat), barrier_stiffness,
        global_mu, default_blend_mu, blend_type);

In Python:

.. code-block:: python

    dhat = 1e-3
    barrier_stiffness = 1.0
    global_mu = 0.6

    collisions = ipctk.NormalCollisions()
    collisions.build(collision_mesh, vertices, dhat)

    # Choose a blend type
    blend_type = ipctk.TangentialCollisions.BlendType.AVG

    # Choose a blend function
    def default_blend_mu(mu0, mu1, blend_type):
        if blend_type == ipctk.TangentialCollisions.BlendType.MIN:
            return min(mu0, mu1)
        elif blend_type == ipctk.TangentialCollisions.BlendType.MAX:
            return max(mu0, mu1)
        elif blend_type == ipctk.TangentialCollisions.BlendType.PRODUCT:
            return mu0 * mu1
        elif blend_type == ipctk.TangentialCollisions.BlendType.HARMONIC_MEAN:
            return 2 * (mu0 * mu1) / (mu0 + mu1)
        elif blend_type == ipctk.TangentialCollisions.BlendType.GEOMETRIC_MEAN:
            return np.sqrt(mu0 * mu1)
        elif blend_type == ipctk.TangentialCollisions.BlendType.AVG:
            return (mu0 + mu1) / 2


    tangential_collisions = ipctk.TangentialCollisions()
    tangential_collisions.build(
        collision_mesh, vertices, collisions, ipctk.BarrierPotential(dhat),
        barrier_stiffness, global_mu, default_blend_mu, blend_type)


Assign Material Properties for Friction
------------------------------------------

Define material IDs and assign friction coefficients (static and kinetic friction values) for pairs of materials. Use `MaterialPairFriction` to specify friction interactions between materials.
The friction coefficients between materials can be found here https://www.engineeringtoolbox.com/friction-coefficients-d_778.html .

.. code-block:: cpp

    global_mu = 0.6;
    global_static_mu = 0.5;
    global_dynamic_mu = 0.3;
    material_pair_friction[{1, 2}] = {0.5, 0.3}; // Friction between material 1 and 2
    material_pair_friction[{2, 3}] = {0.4, 0.25}; // Friction between material 2 and 3


    // Material IDs and corresponding friction parameters
    std::map<std::pair<int, int>, ipc::TangentialCollision::MaterialPairFriction> material_pair_friction;
    material_pair_friction[{1, 2}] = {0.5, 0.3}; // Friction between material 1 and 2
    material_pair_friction[{2, 3}] = {0.4, 0.25}; // Friction between material 2 and 3


    tangential_collisions = ipctk.TangentialCollisions()
    tangential_collisions.build(
        collision_mesh, vertices, collisions, ipctk.BarrierPotential(dhat),
        barrier_stiffness, global_mu, global_static_mu, global_dynamic_mu, material_pair_friction)


.. code-block:: python

    global_mu = 0.6
    global_static_mu = 0.5
    global_dynamic_mu = 0.3
    material_pair_friction = {
        (1, 2): (0.5, 0.3),  # Friction between material 1 and 2
        (2, 3): (0.4, 0.25)  # Friction between material 2 and 3
    }

    tangential_collisions = ipctk.TangentialCollisions()
    tangential_collisions.build(
        collision_mesh, vertices, collisions, ipctk.BarrierPotential(dhat),
        barrier_stiffness, global_mu, global_static_mu, global_dynamic_mu, material_pair_friction)


Compute the Friction Dissipative Potential
---------------------------------------------

Use the ``FrictionPotential`` class to calculate the friction dissipative potential. Define `eps_v` as the threshold for static friction.

.. code-block:: cpp

    const double eps_v = 1e-3;
    ipc::FrictionPotential D(eps_v);

    Eigen::MatrixXd velocity = vertices - collision_mesh.rest_positions();
    double friction_potential = D(tangential_collisions, collision_mesh, velocity);

In Python:

.. code-block:: python

    eps_v = 1e-3
    D = ipctk.FrictionPotential(eps_v)

    velocity = vertices - collision_mesh.rest_positions
    friction_potential = D(tangential_collisions, collision_mesh, velocity)

Compute Friction Gradients and Hessians
------------------------------------------

The gradient and Hessian of the friction dissipative potential can be computed to model frictional forces in iterative solvers.

.. code-block:: cpp

    Eigen::VectorXd friction_potential_grad = D.gradient(tangential_collisions, collision_mesh, velocity);
    Eigen::SparseMatrix<double> friction_potential_hess = D.hessian(tangential_collisions, collision_mesh, velocity);

In Python:

.. code-block:: python

    friction_potential_grad = D.gradient(tangential_collisions, collision_mesh, velocity)
    friction_potential_hess = D.hessian(tangential_collisions, collision_mesh, velocity)