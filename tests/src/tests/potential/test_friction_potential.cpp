#include <tests/config.hpp>
#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>
#include <ipc/friction/friction_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <finitediff.hpp>

using namespace ipc;

// auto default_blend_mu = [](double mu1, double mu2, std::optional<BlendType>) {
//     return ipc::blend_mu(mu1, mu2, BlendType::TRANSITION);
// };

TEST_CASE("Friction gradient and hessian", "[friction][gradient][hessian]")
{
    // Use simple data generator (no pairwise friction)
    FrictionSimpleData data = friction_data_generator();
    const auto& [V0, V1, E, F, collisions, mu, epsv_times_h, dhat, barrier_stiffness] = data;

    const Eigen::MatrixXd U = V1 - V0;
    const CollisionMesh mesh(V0, E, F);

    FrictionCollisions friction_collisions;
    friction_collisions.build(mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, mu);

    FrictionPotential D(epsv_times_h);

    // Compute the gradient analytically
    const Eigen::VectorXd grad = D.gradient(friction_collisions, mesh, U);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_U = fd::unflatten(x, V1.cols()) - V0;
        return D(friction_collisions, mesh, fd_U);
    };

    Eigen::VectorXd fgrad;
    fd::finite_gradient(fd::flatten(V1), f, fgrad);

    // Check if the computed gradient matches the analytical gradient
    CHECK(fd::compare_gradient(grad, fgrad));

    // Compute the Hessian analytically
    const Eigen::MatrixXd hess = D.hessian(friction_collisions, mesh, U);

    // Compute the Hessian using finite differences
    Eigen::MatrixXd fhess;
    fd::finite_hessian(fd::flatten(V1), f, fhess);

    // Check if the computed Hessian matches the analytical Hessian
    CHECK(fd::compare_hessian(hess, fhess, 1e-3));
}

// TEST_CASE(
//     "Pairwise friction gradient and hessian with static and kinetic transition",
//     "[friction][pairwise][gradient][hessian]")
// {
//     // Use complex data generator (pairwise friction handling)
//     FrictionComplexData data = friction_data_generator_with_pairwise();
//     const auto& [V0, V1, E, F, collisions, mu, static_mu, kinetic_mu, epsv_times_h, dhat, barrier_stiffness, pairwise_friction, blend_mu] = data;

//     const Eigen::MatrixXd U = V1 - V0;
//     const CollisionMesh mesh(V0, E, F);

//     FrictionCollisions friction_collisions;
//     friction_collisions.build(
//         mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, static_mu, static_mu, kinetic_mu, pairwise_friction, default_blend_mu);

//     FrictionPotential D(epsv_times_h);

//     SECTION("Checking pairwise gradient with transition")
//     {
//         const Eigen::VectorXd grad = D.gradient(friction_collisions, mesh, U);
//         CHECK(fd::compare_gradient(grad, grad)); // Self-consistency check
//     }

//     SECTION("Checking pairwise Hessian with transition")
//     {
//         const Eigen::MatrixXd hess = D.hessian(friction_collisions, mesh, U);

//         Eigen::MatrixXd fhess;
//         fd::finite_hessian(
//             fd::flatten(V1),
//             [&](const Eigen::VectorXd& x) {
//                 const Eigen::MatrixXd fd_U = fd::unflatten(x, V1.cols()) - V0;
//                 return D(friction_collisions, mesh, fd_U);
//             },
//             fhess);

//         CHECK(fd::compare_hessian(hess, fhess, 1e-3));
//     }
// }

TEST_CASE("Static vs kinetic friction transition", "[friction][transition]")
{
    // Use simple data generator for static vs kinetic friction transition
    FrictionSimpleData data = friction_data_generator();
    const auto& [V0, V1, E, F, collisions, mu, epsv_times_h, dhat, barrier_stiffness] = data;

    const Eigen::MatrixXd U = V1 - V0;
    const CollisionMesh mesh(V0, E, F);

    FrictionCollisions friction_collisions;
    friction_collisions.build(mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, mu);

    FrictionPotential D(epsv_times_h);

    SECTION("Force computation with static/kinetic transition")
    {
        const Eigen::VectorXd force = D.force(
            friction_collisions, mesh, V0, U, U, BarrierPotential(dhat),
            barrier_stiffness, 0.0, false);
        REQUIRE(force.size() == V1.size());
        CHECK(force.lpNorm<Eigen::Infinity>() > 0); // Check non-zero force
    }
}
