#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <ipc/friction/friction_constraints.hpp>
#include <ipc/potentials/friction_potential.hpp>

using namespace ipc;

TEST_CASE("Friction Potential Refactor", "[potential][friction_potential]")
{
    FrictionData data = friction_data_generator();
    const auto& [vertices_t0, vertices_t1, edges, faces, collision_constraints, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;

    const Eigen::MatrixXd velocities = vertices_t1 - vertices_t0;

    const CollisionMesh mesh(vertices_t0, edges, faces);

    FrictionConstraints contacts;
    contacts.build(
        mesh, vertices_t0, collision_constraints, dhat, barrier_stiffness, mu);

    const FrictionPotential D(epsv_times_h);
}
