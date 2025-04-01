#include <tests/config.hpp>
#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/ipc.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/potentials/barrier_potential.hpp>

using namespace ipc;

TEST_CASE("Simulation with material friction", "[friction][material][simulation]")
{
    // Create a simple collision scene with two objects
    // Object 1: Material 0
    Eigen::MatrixXd V1(4, 3);
    V1 << -0.5, -0.5, 0.0,
           0.5, -0.5, 0.0,
           0.5,  0.5, 0.0,
          -0.5,  0.5, 0.0; // Square in XY plane
          
    Eigen::MatrixXi E1(4, 2);
    E1 << 0, 1, 1, 2, 2, 3, 3, 0;
    
    // Object 2: Material 1
    Eigen::MatrixXd V2 = V1;
    V2.col(2).array() += 0.01; // Lift slightly above object 1
    V2.col(0).array() += 0.5;  // Shift to the right
    Eigen::MatrixXi E2 = E1;
    
    // Combine into a single mesh
    Eigen::MatrixXd V(8, 3);
    V << V1, V2;
    
    Eigen::MatrixXi E(8, 2);
    for (int i = 0; i < 4; i++) {
        E.row(i) = E1.row(i);
        E.row(i+4) << E2(i, 0) + 4, E2(i, 1) + 4; // Offset indices for obj 2
    }
    
    // Create collision mesh
    CollisionMesh mesh(V, E);
    
    // Set material IDs
    Eigen::VectorXi vertex_materials(8);
    vertex_materials << 0, 0, 0, 0, 1, 1, 1, 1;  // First 4 vertices are material 0, last 4 are material 1
    mesh.set_vertex_materials(vertex_materials);
    
    // Create a material friction lookup table
    auto friction_table = std::make_shared<Eigen::MatrixXd>(2, 2);
    (*friction_table) << 
        0.2, 0.5,  // material 0 against 0 or 1
        0.5, 0.3;  // material 1 against 0 or 1
    
    // Create a friction potential
    FrictionPotential friction_potential(1e-5);
    friction_potential.set_material_friction_table(friction_table);
    
    // Apply a small displacement to simulate objects moving
    Eigen::MatrixXd displacement = Eigen::MatrixXd::Zero(V.rows(), 3);
    displacement.block(4, 0, 4, 2) = Eigen::MatrixXd::Ones(4, 2) * 0.01; // Move object 2
    
    // Create normal collisions
    NormalCollisions normal_collisions;
    normal_collisions.set_enable_shape_derivatives(true);
    normal_collisions.build(mesh, V + displacement, 0.1);
    
    // Create tangential collisions
    TangentialCollisions tangential_collisions;
    tangential_collisions.build(mesh, V + displacement, normal_collisions, 
                                BarrierPotential(0.1), 100.0);
    
    // Verify the collisions were created
    REQUIRE(!normal_collisions.empty());
    REQUIRE(!tangential_collisions.empty());
    
    // Check that material IDs were properly assigned
    for (size_t i = 0; i < tangential_collisions.size(); i++) {
        const TangentialCollision& collision = tangential_collisions[i];
        CHECK(collision.has_material_ids());
        
        // Validate material IDs make sense for this collision
        int mat1 = collision.material_id1;
        int mat2 = collision.material_id2;
        
        // Material IDs should be either 0 or 1
        CHECK((mat1 == 0 || mat1 == 1));
        CHECK((mat2 == 0 || mat2 == 1));
    }
    
    // Apply velocities to test friction
    Eigen::MatrixXd velocities = Eigen::MatrixXd::Zero(V.rows(), 3);
    velocities.block(4, 0, 4, 1) = Eigen::VectorXd::Ones(4) * 0.1; // X velocity for object 2
    
    // Compute friction forces
    Eigen::VectorXd friction_force = friction_potential.force(
        tangential_collisions, mesh, V, displacement, velocities,
        BarrierPotential(0.1), 100.0, 0, false);
    
    // The force should be non-zero
    CHECK(friction_force.norm() > 0);
    
    // Change the friction table to see if forces change
    (*friction_table) << 
        0.1, 0.1,  // material 0 against 0 or 1 (reduced friction)
        0.1, 0.1;  // material 1 against 0 or 1 (reduced friction)
    
    Eigen::VectorXd reduced_friction_force = friction_potential.force(
        tangential_collisions, mesh, V, displacement, velocities,
        BarrierPotential(0.1), 100.0, 0, false);
    
    // The force should be less than before
    CHECK(reduced_friction_force.norm() < friction_force.norm());
}
