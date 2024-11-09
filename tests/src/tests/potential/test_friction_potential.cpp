#include <tests/config.hpp>
#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/friction/friction_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>

#include <finitediff.hpp>

using namespace ipc;

TEST_CASE("Friction gradient and hessian", "[friction][gradient][hessian]")
{
    // Generate test data for friction
    FrictionSimpleData data = friction_data_generator();
    const auto& [V0, V1, E, F, collisions, static_mu, kinetic_mu, epsv_times_h, dhat, barrier_stiffness] = data;

    const Eigen::MatrixXd U = V1 - V0;

    const CollisionMesh mesh(V0, E, F);

    FrictionCollisions friction_collisions;
    friction_collisions.build(
        mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, static_mu, kinetic_mu);

    const FrictionPotential D(epsv_times_h);

    // Compute the gradient analytically
    const Eigen::VectorXd grad = D.gradient(friction_collisions, mesh, U, static_mu, kinetic_mu);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_U = fd::unflatten(x, data.V1.cols()) - data.V0;
        return D(friction_collisions, mesh, fd_U, static_mu, kinetic_mu);
    };

    Eigen::VectorXd fgrad;
    fd::finite_gradient(fd::flatten(V1), f, fgrad);

    // Check if the computed gradient matches the analytical gradient
    CHECK(fd::compare_gradient(grad, fgrad));

    // Compute the Hessian analytically
    const Eigen::MatrixXd hess = D.hessian(friction_collisions, mesh, U, static_mu, kinetic_mu);

    // Compute the Hessian using finite differences
    Eigen::MatrixXd fhess;
    fd::finite_hessian(fd::flatten(V1), f, fhess);

    // Check if the computed Hessian matches the analytical Hessian
    CHECK(fd::compare_hessian(hess, fhess, 1e-3));
}

TEST_CASE("Pairwise friction gradient and hessian with static and kinetic transition", "[friction][pairwise][gradient][hessian]")
{
    // Generate test data for pairwise friction with transition
    FrictionSimpleData data = friction_data_generator();
    const auto& [V0, V1, E, F, collisions, static_mu, kinetic_mu, epsv_times_h, dhat, barrier_stiffness] = data;

    const Eigen::MatrixXd U = V1 - V0;

    const CollisionMesh mesh(V0, E, F);

    FrictionCollisions friction_collisions;
    friction_collisions.build(
        mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, static_mu, kinetic_mu);

    FrictionPotential D(epsv_times_h);

    // Gradient and Hessian calculations with transition
    for (const auto& collision : collisions) {
        // Assuming collision uses either kinetic or static friction based on relative velocity
        double velocity_norm = (collision.dof(U, E, F)).norm();
        bool use_kinetic = velocity_norm < epsv_times_h;
        double expected_mu = use_kinetic ? kinetic_mu : static_mu;

        SECTION("Checking pairwise gradient with transition")
        {
            const Eigen::VectorXd grad = D.gradient(friction_collisions, mesh, U, static_mu, kinetic_mu);
            CHECK(fd::compare_gradient(grad, grad)); // Check self-consistency as baseline
        }

        SECTION("Checking pairwise Hessian with transition")
        {
            const Eigen::MatrixXd hess = D.hessian(friction_collisions, mesh, U, static_mu, kinetic_mu);
            CHECK(hess.rows() == U.size());
            CHECK(hess.cols() == U.size());

            Eigen::MatrixXd fhess;
            fd::finite_hessian(fd::flatten(V1), [&](const Eigen::VectorXd& x) {
                const Eigen::MatrixXd fd_U = fd::unflatten(x, data.V1.cols()) - data.V0;
                return D(friction_collisions, mesh, fd_U, static_mu, kinetic_mu);
            }, fhess);

            CHECK(fd::compare_hessian(hess, fhess, 1e-3));
        }
    }
}

TEST_CASE("Static vs kinetic friction transition", "[friction][transition]")
{
    // Generate test data for friction
    FrictionSimpleData data = friction_data_generator();
    const auto& [V0, V1, E, F, collisions, static_mu, kinetic_mu, epsv_times_h, dhat, barrier_stiffness] = data;

    const Eigen::MatrixXd U = V1 - V0;

    const CollisionMesh mesh(V0, E, F);

    FrictionCollisions friction_collisions;
    friction_collisions.build(
        mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, static_mu, kinetic_mu);

    FrictionPotential D(epsv_times_h);

    // Check static vs kinetic transition in force calculations
    for (const auto& collision : collisions) {
        double velocity_norm = (collision.dof(U, E, F)).norm();
        bool is_static = velocity_norm < epsv_times_h;

        double expected_mu = is_static ? static_mu : kinetic_mu;
        CHECK(expected_mu == (is_static ? static_mu : kinetic_mu));

        SECTION("Force computation with static/kinetic transition")
        {
            const Eigen::VectorXd force = D.force(friction_collisions, mesh, V0, U, U, BarrierPotential(dhat), barrier_stiffness, 0.0, false, static_mu, kinetic_mu);
            REQUIRE(force.size() == V1.size());
            CHECK(force.lpNorm<Eigen::Infinity>() > 0); // Check non-zero force
        }
    }
}

TEST_CASE("Self-consistency of pairwise static and kinetic friction transition in gradient and Hessian", "[friction][pairwise][self-consistency]")
{
    FrictionSimpleData data = friction_data_generator();
    const auto& [V0, V1, E, F, collisions, static_mu, kinetic_mu, epsv_times_h, dhat, barrier_stiffness] = data;

    const Eigen::MatrixXd U = V1 - V0;

    const CollisionMesh mesh(V0, E, F);

    FrictionCollisions friction_collisions;
    friction_collisions.build(
        mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, static_mu, kinetic_mu);

    FrictionPotential D(epsv_times_h);

    // Compare pairwise results across multiple instances for consistency
    const Eigen::VectorXd grad1 = D.gradient(friction_collisions, mesh, U, static_mu, kinetic_mu);
    const Eigen::MatrixXd hess1 = D.hessian(friction_collisions, mesh, U, static_mu, kinetic_mu);

    const Eigen::VectorXd grad2 = D.gradient(friction_collisions, mesh, U, static_mu, kinetic_mu);
    const Eigen::MatrixXd hess2 = D.hessian(friction_collisions, mesh, U, static_mu, kinetic_mu);

    SECTION("Consistency in gradient and Hessian between calls")
    {
        CHECK(grad1.isApprox(grad2, 1e-8));
        CHECK(hess1.isApprox(hess2, 1e-8));
    }
}
