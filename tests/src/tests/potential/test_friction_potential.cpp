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
    const auto& [V0, V1, E, F, collisions, mu, epsv_times_h, dhat, barrier_stiffness] = data;

    const Eigen::MatrixXd U = V1 - V0;

    const CollisionMesh mesh(V0, E, F);

    FrictionCollisions friction_collisions;
    friction_collisions.build(mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, mu);

    const FrictionPotential D(epsv_times_h);

    // Compute the gradient analytically
    const Eigen::VectorXd grad = D.gradient(friction_collisions, mesh, U);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_U = fd::unflatten(x, data.V1.cols()) - data.V0;
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

// TEST_CASE("Pairwise friction hessian", "[friction][pairwise][hessian]")
// {
//     // Generate complex test data for pairwise friction
//     FrictionComplexData data = friction_data_generator_with_pairwise();
//     const auto& [V0, V1, E, F, collisions, static_mu, kinetic_mu, epsv_times_h, dhat, barrier_stiffness, pairwise_friction] = data;

//     const Eigen::MatrixXd U = V1 - V0;

//     const CollisionMesh mesh(V0, E, F);

//     FrictionCollisions friction_collisions;
//     friction_collisions.build(mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, static_mu, kinetic_mu, pairwise_friction);

//     const FrictionPotential D(epsv_times_h);

//     // Compute the pairwise friction Hessian using the analytical method
//     const MatrixMax12d hess = D.pairwise_hessian(friction_collisions, V1, FrictionPotential::PSDProjectionMethod::NONE, static_mu, kinetic_mu);

//     // Check Hessian with finite differences
//     auto f = [&](const Eigen::VectorXd& x) {
//         const Eigen::MatrixXd fd_U = fd::unflatten(x, data.V1.cols()) - data.V0;
//         return D.pairwise_hessian(friction_collisions, fd_U, FrictionPotential::PSDProjectionMethod::NONE, static_mu, kinetic_mu);
//     };

//     Eigen::MatrixXd fd_hessian;
//     fd::finite_hessian(fd::flatten(V1), f, fd_hessian);

//     // Check if the computed pairwise Hessian matches the analytical Hessian
//     CHECK(fd::compare_hessian(hess, fd_hessian, 1e-3));
// }
