#include <tests/config.hpp>
#include <tests/friction/friction_data_generator.hpp>
#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/collisions/tangential/tangential_collisions.hpp>
#include <ipc/potentials/friction_potential.hpp>
#include <ipc/potentials/barrier_potential.hpp>


#include <finitediff.hpp>

using namespace ipc;

TEST_CASE("Friction gradient and hessian", "[friction][gradient][hessian]")
{
    FrictionData data = friction_data_generator();
    const auto& [V0, V1, E, F, collisions, mu, epsv_times_h, dhat, barrier_stiffness] =
        data;

    const Eigen::MatrixXd U = V1 - V0;

    const CollisionMesh mesh(V0, E, F);

    TangentialCollisions tangential_collisions;
    tangential_collisions.build(
        mesh, V0, collisions, BarrierPotential(dhat), barrier_stiffness, mu);

    const FrictionPotential D(epsv_times_h);

    const Eigen::VectorXd grad = D.gradient(tangential_collisions, mesh, U);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        const Eigen::MatrixXd fd_U = fd::unflatten(x, data.V1.cols()) - data.V0;
        return D(tangential_collisions, mesh, fd_U);
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(fd::flatten(V1), f, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));

    const Eigen::MatrixXd hess = D.hessian(tangential_collisions, mesh, U);
    Eigen::MatrixXd fhess;
    fd::finite_hessian(fd::flatten(V1), f, fhess);
    CHECK(fd::compare_hessian(hess, fhess, 1e-3));
}