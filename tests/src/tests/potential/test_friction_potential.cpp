#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/friction/friction_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>

using namespace ipc;

TEST_CASE("Friction Potential Refactor", "[potential][friction_potential]")
{
    FrictionData data = friction_data_generator();
    const auto& [vertices_t0, vertices_t1, edges, faces, collisions, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;

    const Eigen::MatrixXd rest_positions =
        vertices_t0 - (vertices_t1 - vertices_t0);
    const Eigen::MatrixXd lagged_displacements = vertices_t0 - rest_positions;
    const Eigen::MatrixXd velocities = vertices_t1 - vertices_t0;

    const CollisionMesh mesh(rest_positions, edges, faces);

    if (collisions.compute_minimum_distance(mesh, rest_positions) == 0
        || collisions.compute_minimum_distance(mesh, vertices_t0) == 0
        || collisions.compute_minimum_distance(mesh, vertices_t1) == 0) {
        return;
    }

    FrictionCollisions friction_collisions;
    friction_collisions.build(
        mesh, vertices_t0, collisions, dhat, barrier_stiffness, mu);

    const FrictionPotential D(epsv_times_h);
}
