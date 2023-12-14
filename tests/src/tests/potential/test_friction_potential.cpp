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

    FrictionPotential D(epsv_times_h);

    const double expected_potential =
        contacts.compute_potential(mesh, velocities, epsv_times_h);
    const double actual_potential = D(mesh, velocities, contacts);
    CHECK(actual_potential == Catch::Approx(expected_potential));

    const Eigen::VectorXd expected_gradient =
        contacts.compute_potential_gradient(mesh, velocities, epsv_times_h);
    const Eigen::VectorXd actual_gradient =
        D.gradient(mesh, velocities, contacts);
    CHECK(actual_gradient.isApprox(expected_gradient));

    Eigen::SparseMatrix<double> expected_hessian, actual_hessian;
    expected_hessian = contacts.compute_potential_hessian(
        mesh, velocities, epsv_times_h, false);
    actual_hessian = D.hessian(mesh, velocities, contacts, false);
    CHECK(actual_hessian.isApprox(expected_hessian));

    // Projected hessian
    expected_hessian = contacts.compute_potential_hessian(
        mesh, velocities, epsv_times_h, true);
    actual_hessian = D.hessian(mesh, velocities, contacts, true);
    INFO("project_hessian_to_psd = true");
    CHECK(actual_hessian.isApprox(expected_hessian));
}
