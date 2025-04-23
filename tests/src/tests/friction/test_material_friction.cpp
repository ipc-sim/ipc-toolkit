#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/ipc.hpp>
#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/potentials/barrier_potential.hpp>

using namespace ipc;

TEST_CASE("Material-specific friction", "[friction][material]")
{
    // Create a simple scene with two objects of different materials
    Eigen::MatrixXd V(4, 3);
    V << 0, 0, 0, 
         1, 0, 0, 
         0, 1, 0,
         1, 1, 0; // 2D square in XY plane
         
    Eigen::MatrixXi E(4, 2);
    E << 0, 1,
         1, 3,
         3, 2,
         2, 0; // square edges
         
    CollisionMesh mesh(V, E);
    
    // Set material IDs - half of the vertices are material 0, half material 1
    Eigen::VectorXi vertex_materials(4);
    vertex_materials << 0, 0, 1, 1;
    mesh.set_vertex_materials(vertex_materials);
    CHECK(mesh.has_material_ids());
    
    // Edge materials derived from vertices
    Eigen::VectorXi edge_materials(4);
    edge_materials << 0, 0, 1, 1;
    mesh.set_edge_materials(edge_materials);
    
    // Simulate a vertex-vertex collision between material 0 and 1
    VertexVertexNormalCollision normal_collision(0, 2);
    normal_collision.material_id1 = mesh.vertex_material(0); // should be 0
    normal_collision.material_id2 = mesh.vertex_material(2); // should be 1
    CHECK(normal_collision.material_id1 == 0);
    CHECK(normal_collision.material_id2 == 1);
    
    // Create a friction collision from the normal collision
    VertexVertexTangentialCollision friction_collision(normal_collision);
    CHECK(friction_collision.material_id1 == 0);
    CHECK(friction_collision.material_id2 == 1);
    CHECK(friction_collision.has_material_ids());
    
    // Create a material-specific friction table
    Eigen::MatrixXd friction_table(2, 2);
    friction_table << 0.3, 0.5,
                      0.5, 0.7; // Different mu values for different material pairs
    
    // Set default friction coefficient
    friction_collision.mu = 0.4;
    
    // Create friction potential
    FrictionPotential friction_potential(0.001);
    
    // Test friction coefficient lookup
    double mu = friction_potential.get_friction_coefficient(friction_collision, friction_table);
    CHECK(mu == Catch::Approx(0.5)); // Should be the value for materials 0,1
    
    // Test with a different material combination
    friction_collision.material_id1 = 1;
    friction_collision.material_id2 = 1;
    mu = friction_potential.get_friction_coefficient(friction_collision, friction_table);
    CHECK(mu == Catch::Approx(0.7)); // Should be the value for materials 1,1
    
    // Test with invalid material IDs
    friction_collision.material_id1 = -1; // NO_MATERIAL_ID
    mu = friction_potential.get_friction_coefficient(friction_collision, friction_table);
    CHECK(mu == Catch::Approx(friction_collision.mu)); // Should fall back to collision's mu
    
    // Set valid material IDs but try with empty friction table
    friction_collision.material_id1 = 0;
    friction_collision.material_id2 = 1;
    mu = friction_potential.get_friction_coefficient(friction_collision, Eigen::MatrixXd(0, 0));
    CHECK(mu == Catch::Approx(friction_collision.mu)); // Should fall back to collision's mu
    
    // Test with material IDs out of range
    friction_collision.material_id1 = 5; // Out of range
    mu = friction_potential.get_friction_coefficient(friction_collision, friction_table);
    CHECK(mu == Catch::Approx(friction_collision.mu)); // Should fall back to collision's mu
}
