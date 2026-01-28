#include <tests/utils.hpp>

#include <catch2/catch_test_macros.hpp>

#include <ipc/dynamics/affine/orthogonality_potential.hpp>

#include <finitediff.hpp>

TEST_CASE("Orthogonality potential", "[abd]")
{
    using namespace ipc;
    using namespace ipc::affine;

    AffineBody body(MatrixMax9d::Random(3, 3), VectorMax9d::Random(3), 1);

    VectorMax12d dof(12);
    dof.head(3) = body.p;
    dof.segment(3, 3) = body.A.col(0);
    dof.segment(6, 3) = body.A.col(1);
    dof.segment(9, 3) = body.A.col(2);

    OrthogonalityPotential V_perp(1);

    // -------------------------------------------------------------------------
    // Gradient
    // -------------------------------------------------------------------------

    const Eigen::VectorXd grad_V_perp = V_perp.gradient(body);
    REQUIRE(grad_V_perp.squaredNorm() > 1e-8);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad_V_perp;
    {
        auto f = [&V_perp, &body](const Eigen::VectorXd& x) {
            AffineBody fd_body;
            fd_body.p = x.head(3);
            fd_body.A.resize(3, 3);
            fd_body.A.col(0) = x.segment(3, 3);
            fd_body.A.col(1) = x.segment(6, 3);
            fd_body.A.col(2) = x.segment(9, 3);
            fd_body.volume = body.volume;
            return V_perp(fd_body);
        };
        fd::finite_gradient(dof, f, fgrad_V_perp);
    }

    CHECK(fd::compare_gradient(grad_V_perp, fgrad_V_perp));
    if (!fd::compare_gradient(grad_V_perp, fgrad_V_perp)) {
        logger().error("grad_V_perp:  {}", grad_V_perp);
        logger().error("fgrad_V_perp: {}", fgrad_V_perp);
    }

    // -------------------------------------------------------------------------
    // Hessian
    // -------------------------------------------------------------------------

    const Eigen::MatrixXd hess_V_perp = V_perp.hessian(body);
    REQUIRE(hess_V_perp.squaredNorm() > 1e-3);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fhess_V_perp;
    {
        auto f = [&](const Eigen::VectorXd& x) {
            AffineBody fd_body;
            fd_body.p = x.head(3);
            fd_body.A.resize(3, 3);
            fd_body.A.col(0) = x.segment(3, 3);
            fd_body.A.col(1) = x.segment(6, 3);
            fd_body.A.col(2) = x.segment(9, 3);
            fd_body.volume = body.volume;
            return V_perp.gradient(fd_body);
        };
        fd::finite_jacobian(dof, f, fhess_V_perp);
    }

    CHECK(fd::compare_hessian(hess_V_perp, fhess_V_perp, 1e-3));
    if (!fd::compare_hessian(hess_V_perp, fhess_V_perp, 1e-3)) {
        tests::print_compare_nonzero(hess_V_perp, fhess_V_perp);
    }
}