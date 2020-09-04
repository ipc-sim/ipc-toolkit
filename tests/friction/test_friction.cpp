#include <catch2/catch.hpp>

#include <finitediff.hpp>
#include <igl/edges.h>

#include <friction/friction.hpp>

#include "../test_utils.hpp"

using namespace ipc;

TEST_CASE("Test friction gradient", "[friction][grad][hess]")
{
    double mu = GENERATE(range(0.0, 1.0, 0.1));
    double epsv_times_h_squared = pow(10, -2 * GENERATE(range(0, 6)));
    double dhat_squared = pow(10, -2 * GENERATE(range(0, 4)));
    double dhat = sqrt(dhat_squared);
    double barrier_stiffness = pow(10, GENERATE(range(0, 2)));
    CAPTURE(mu, epsv_times_h_squared, dhat_squared, barrier_stiffness);

    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi E, F;
    Candidates contact_constraint_set;
    SECTION("point-triangle")
    {
        V0.resize(4, 3);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << 0, d, 0;   // point at t=0
        V0.row(1) << -1, 0, 1;  // triangle vertex 0 at t=0
        V0.row(2) << 2, 0, 0;   // triangle vertex 1 at t=0
        V0.row(3) << -1, 0, -1; // triangle vertex 2 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V1.row(0) << 1, d + dy, 0; // point at t=1

        F.resize(1, 3);
        F << 1, 2, 3;
        igl::edges(F, E);
        REQUIRE(E.rows() == 3);

        contact_constraint_set.fv_candidates.emplace_back(0, 0);
    }
    SECTION("edge-edge")
    {
        V0.resize(4, 3);
        double d = GENERATE_COPY(range(0.0, 2 * dhat, 2 * dhat / 10.0));
        V0.row(0) << -1, d, 0; // edge a vertex 0 at t=0
        V0.row(1) << 1, d, 0;  // edge a vertex 1 at t=0
        V0.row(2) << 0, 0, -1; // edge b vertex 0 at t=0
        V0.row(3) << 0, 0, 1;  // edge b vertex 1 at t=0

        V1 = V0;
        double dy = GENERATE(-1, 1, 1e-1);
        V0.row(0) << 0.5, d, 0; // edge a vertex 0 at t=1
        V0.row(1) << 2.5, d, 0; // edge a vertex 1 at t=1

        E.resize(2, 2);
        E.row(0) << 0, 1;
        E.row(1) << 2, 3;

        contact_constraint_set.ee_candidates.emplace_back(0, 1);
    }

    Candidates friction_constraint_set;
    std::vector<Eigen::VectorXd> closest_points;
    std::vector<Eigen::MatrixXd> tangent_bases;
    Eigen::VectorXd normal_force_magnitudes;
    compute_friction_bases(
        V0, E, F, contact_constraint_set, dhat_squared, barrier_stiffness,
        friction_constraint_set, closest_points, tangent_bases,
        normal_force_magnitudes);

    Eigen::VectorXd grad = compute_friction_potential_gradient(
        V0, V1, E, F, friction_constraint_set, closest_points, tangent_bases,
        normal_force_magnitudes, epsv_times_h_squared, mu);

    // Compute the gradient using finite differences
    auto f = [&](const Eigen::VectorXd& x) {
        return ipc::compute_friction_potential(
            V0, unflatten(x, V1.cols()), E, F, friction_constraint_set,
            closest_points, tangent_bases, normal_force_magnitudes,
            epsv_times_h_squared, mu);
    };
    Eigen::VectorXd fgrad;
    fd::finite_gradient(flatten(V1), f, fgrad);
    CHECK(fd::compare_gradient(grad, fgrad));

    Eigen::MatrixXd hess = compute_friction_potential_hessian(
        V0, V1, E, F, friction_constraint_set, closest_points, tangent_bases,
        normal_force_magnitudes, epsv_times_h_squared, mu);
    Eigen::MatrixXd fhess;
    fd::finite_hessian(flatten(V1), f, fhess);
    if (!fd::compare_hessian(hess, fhess, 1e-3)) {
        std::cout << hess.transpose() << std::endl;
        std::cout << fhess.transpose() << std::endl;
    }
    CHECK(fd::compare_hessian(hess, fhess, 1e-3));
}
