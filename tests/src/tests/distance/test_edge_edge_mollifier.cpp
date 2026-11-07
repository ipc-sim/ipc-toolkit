#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>
#include <ipc/potentials/barrier_potential.hpp>

#include <finitediff.hpp>

using namespace ipc;

std::array<Eigen::Vector3d, 4> get_edges()
{
    Eigen::Vector3d ea0, ea1, eb0, eb1;
    SECTION("Perp")
    {
        ea0 = Eigen::Vector3d(-1, 0, 0);
        ea1 = Eigen::Vector3d(1, 0, 0);
        eb0 = Eigen::Vector3d(0, -1, 0);
        eb1 = Eigen::Vector3d(0, 1, 0);
        CHECK(
            edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1)
            == Catch::Approx(16));
    }
    SECTION("Almost parallel")
    {
        ea0 = Eigen::Vector3d(-1, 0, 0);
        ea1 = Eigen::Vector3d(1, 0, 0);
        eb0 = Eigen::Vector3d(-1, 1e-9, 0);
        eb1 = Eigen::Vector3d(1, -1e-9, 0);
        CHECK(
            edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1)
            == Catch::Approx(0).margin(1e-9));
    }
    return { { ea0, ea1, eb0, eb1 } };
}

TEST_CASE("Edge-Edge Cross Squarednorm", "[distance][edge-edge][mollifier]")
{
    const auto& [ea0, ea1, eb0, eb1] = get_edges();

    Vector12d x;
    x << ea0, ea1, eb0, eb1;

    const Vector12d grad =
        edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& _x) {
            return edge_edge_cross_squarednorm(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                _x.segment<3>(9));
        },
        fgrad);

    CHECK(fd::compare_gradient(grad, fgrad));

    const Matrix12d hess =
        edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1);

    // Compute the hessian using finite differences
    Eigen::MatrixXd fhess;
    fd::finite_hessian(
        x,
        [](const Eigen::VectorXd& _x) {
            return edge_edge_cross_squarednorm(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                _x.segment<3>(9));
        },
        fhess);

    CHECK(fd::compare_hessian(hess, fhess));
}

TEST_CASE("Edge-Edge Mollifier Scalar", "[distance][edge-edge][mollifier]")
{
    const double rel_x = GENERATE(0, 0.5, 1, 2);
    const double eps_x = GENERATE(1e-3, 1e-1, 1, 2);
    const double x = rel_x * eps_x;
    CAPTURE(rel_x, eps_x, x);

    const double m = edge_edge_mollifier(x, eps_x);
    CHECK(m >= 0);
    CHECK(m <= 1);
    if (x > eps_x)
        CHECK(m == 1);

    // ---

    const double dm_dx = edge_edge_mollifier_gradient(x, eps_x);

    // Compute the gradient using finite differences
    Eigen::VectorXd f_dm_dx;
    fd::finite_gradient(
        Vector1d::Constant(x),
        [eps_x](const Eigen::VectorXd& _x) {
            return edge_edge_mollifier(_x(0), eps_x);
        },
        f_dm_dx, fd::SECOND, 1e-12);

    CHECK(dm_dx == Catch::Approx(f_dm_dx(0)));

    // ---

    if (std::abs(x - eps_x) > 1e-6) { // mollifier is only C1 at x = eps_x
        const double d2m_dx2 = edge_edge_mollifier_hessian(x, eps_x);

        // Compute the gradient using finite differences
        Eigen::MatrixXd f_d2m_dx2;
        fd::finite_hessian(
            Vector1d::Constant(x),
            [eps_x](const Eigen::VectorXd& _x) {
                return edge_edge_mollifier(_x(0), eps_x);
            },
            f_d2m_dx2);

        CHECK(d2m_dx2 == Catch::Approx(f_d2m_dx2(0)));
    }

    // ---

    const double dm_deps_x = edge_edge_mollifier_derivative_wrt_eps_x(x, eps_x);

    // Compute the gradient using finite differences
    Eigen::VectorXd f_dm_deps_x;
    fd::finite_gradient(
        Vector1d::Constant(eps_x),
        [x](const Eigen::VectorXd& _eps_x) {
            return edge_edge_mollifier(x, _eps_x(0));
        },
        f_dm_deps_x, fd::SECOND, 1e-12);

    CHECK(dm_deps_x == Catch::Approx(f_dm_deps_x(0)));

    // ---

    if (std::abs(x - eps_x) > 1e-6) { // mollifier is only C1 at x = eps_x
        const double d2m_deps_xdx =
            edge_edge_mollifier_gradient_derivative_wrt_eps_x(x, eps_x);

        // Compute the gradient using finite differences
        Eigen::VectorXd f_d2m_deps_xdx;
        fd::finite_gradient(
            Vector1d::Constant(eps_x),
            [x](const Eigen::VectorXd& _eps_x) {
                return edge_edge_mollifier_gradient(x, _eps_x(0));
            },
            f_d2m_deps_xdx);

        CHECK(d2m_deps_xdx == Catch::Approx(f_d2m_deps_xdx(0)));
    }
}

TEST_CASE("Edge-Edge Mollifier", "[distance][edge-edge][mollifier]")
{
    const Eigen::Vector3d rest_ea0(0, 0, 0);
    const Eigen::Vector3d rest_ea1(1, 0, 0);
    const Eigen::Vector3d rest_eb0(0, 0, 0);
    const Eigen::Vector3d rest_eb1(0, 1, 0);
    const auto& [ea0, ea1, eb0, eb1] = get_edges();

    const double eps_x =
        edge_edge_mollifier_threshold(rest_ea0, rest_ea1, rest_eb0, rest_eb1);

    const double m = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    CHECK(m >= 0);
    CHECK(m <= 1);

    Vector12d x;
    x << rest_ea0, rest_ea1, rest_eb0, rest_eb1;

    Vector12d v;
    v << ea0, ea1, eb0, eb1;

    const Vector12d u = v - x;

    // ---

    const Vector12d grad =
        edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        v,
        [eps_x](const Eigen::VectorXd& _v) {
            return edge_edge_mollifier(
                _v.segment<3>(0), _v.segment<3>(3), _v.segment<3>(6),
                _v.segment<3>(9), eps_x);
        },
        fgrad);

    CHECK(fd::compare_gradient(grad, fgrad));

    // ---

    const Matrix12d hess =
        edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x);

    // Compute the hessian using finite differences
    Eigen::MatrixXd fhess;
    fd::finite_hessian(
        v,
        [eps_x](const Eigen::VectorXd& _v) {
            return edge_edge_mollifier(
                _v.segment<3>(0), _v.segment<3>(3), _v.segment<3>(6),
                _v.segment<3>(9), eps_x);
        },
        fhess);

    CHECK(fd::compare_hessian(hess, fhess));

    // ---

    const Vector12d grad_wrt_x = edge_edge_mollifier_gradient_wrt_x(
        rest_ea0, rest_ea1, rest_eb0, rest_eb1, ea0, ea1, eb0, eb1);

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad_wrt_x;
    fd::finite_gradient(
        x,
        [u](const Eigen::VectorXd& _x) {
            const Eigen::VectorXd _v = _x + u;
            return edge_edge_mollifier(
                _v.segment<3>(0), _v.segment<3>(3), _v.segment<3>(6),
                _v.segment<3>(9),
                edge_edge_mollifier_threshold(
                    _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                    _x.segment<3>(9)));
        },
        fgrad_wrt_x);

    CHECK(fd::compare_gradient(grad_wrt_x, fgrad_wrt_x));

    // ---

    const Matrix12d jac_wrt_x_u = edge_edge_mollifier_gradient_jacobian_wrt_x(
        rest_ea0, rest_ea1, rest_eb0, rest_eb1, ea0, ea1, eb0, eb1);

    // Compute the gradient using finite differences
    Eigen::MatrixXd fjac_wrt_x_u;
    fd::finite_jacobian(
        x,
        [u](const Eigen::VectorXd& _x) {
            const Eigen::VectorXd _v = _x + u;
            return edge_edge_mollifier_gradient(
                _v.segment<3>(0), _v.segment<3>(3), _v.segment<3>(6),
                _v.segment<3>(9),
                edge_edge_mollifier_threshold(
                    _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                    _x.segment<3>(9)));
        },
        fjac_wrt_x_u);

    CAPTURE(jac_wrt_x_u, fjac_wrt_x_u);
    CHECK(fd::compare_jacobian(jac_wrt_x_u, fjac_wrt_x_u));
}

TEST_CASE("Edge-Edge Mollifier Threshold", "[distance][edge-edge][mollifier]")
{
    const auto& [ea0, ea1, eb0, eb1] = get_edges();

    const double eps_x = edge_edge_mollifier_threshold(ea0, ea1, eb0, eb1);
    CHECK(eps_x == Catch::Approx(0.016));

    const Vector12d grad =
        edge_edge_mollifier_threshold_gradient(ea0, ea1, eb0, eb1);

    Vector12d x;
    x << ea0, ea1, eb0, eb1;

    // Compute the gradient using finite differences
    Eigen::VectorXd fgrad;
    fd::finite_gradient(
        x,
        [](const Eigen::VectorXd& _x) {
            return edge_edge_mollifier_threshold(
                _x.segment<3>(0), _x.segment<3>(3), _x.segment<3>(6),
                _x.segment<3>(9));
        },
        fgrad);

    CHECK(fd::compare_gradient(grad, fgrad));
}

// Single mollified EE pair, FD-check shape_derivative against gradient.
TEST_CASE(
    "Edge-Edge Mollifier Shape Derivative",
    "[distance][edge-edge][mollifier][shape_derivative]")
{
    const double dhat = 1e-3;
    const bool aw = GENERATE(false, true);
    const double y_sep = GENERATE(9.99e-4, 9.9e-4, 9.5e-4);
    const double z_twist = GENERATE(5e-3, 1e-2, 2e-2);

    Eigen::MatrixXd rest(4, 3);
    rest << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, y_sep, 0.0, 1.0, y_sep, z_twist;
    const Eigen::MatrixXd displaced = rest;
    Eigen::MatrixXi edges(2, 2);
    edges << 0, 1, 2, 3;
    Eigen::MatrixXi faces(0, 3);

    CollisionMesh mesh(rest, edges, faces);
    if (!mesh.are_area_jacobians_initialized())
        mesh.init_area_jacobians();

    NormalCollisions collisions;
    collisions.set_use_area_weighting(aw);
    collisions.set_enable_shape_derivatives(true);
    collisions.build(mesh, displaced, dhat);
    REQUIRE(collisions.size() == 1);
    REQUIRE(collisions.is_edge_edge(0));
    REQUIRE(collisions[0].is_mollified());

    BarrierPotential bp(
        dhat, /*stiffness=*/1e8, /*use_physical_barrier=*/false);

    Eigen::SparseMatrix<double> JF_an_sparse =
        bp.shape_derivative(collisions, mesh, displaced);
    Eigen::MatrixXd JF_an(JF_an_sparse);

    // Column-wise central FD on rest DOFs, fixed collision set.
    const int n_verts = static_cast<int>(rest.rows());
    const int dim = static_cast<int>(rest.cols());
    const int ndof = n_verts * dim;
    const double eps = 1e-7;
    Eigen::MatrixXd JF_fd(ndof, ndof);
    JF_fd.setZero();
    for (int v = 0; v < n_verts; ++v) {
        for (int d = 0; d < dim; ++d) {
            Eigen::MatrixXd rest_plus = rest;
            Eigen::MatrixXd rest_minus = rest;
            rest_plus(v, d) += eps;
            rest_minus(v, d) -= eps;
            const Eigen::MatrixXd disp_plus = rest_plus + (displaced - rest);
            const Eigen::MatrixXd disp_minus = rest_minus + (displaced - rest);
            CollisionMesh mesh_plus(rest_plus, edges, faces);
            if (mesh.are_area_jacobians_initialized())
                mesh_plus.init_area_jacobians();
            CollisionMesh mesh_minus(rest_minus, edges, faces);
            if (mesh.are_area_jacobians_initialized())
                mesh_minus.init_area_jacobians();
            const Eigen::VectorXd F_plus =
                bp.gradient(collisions, mesh_plus, disp_plus);
            const Eigen::VectorXd F_minus =
                bp.gradient(collisions, mesh_minus, disp_minus);
            JF_fd.col(v * dim + d) = (F_plus - F_minus) / (2.0 * eps);
        }
    }

    const double an = JF_an.norm();
    const double fn = JF_fd.norm();
    const double diff = (JF_an - JF_fd).norm();
    const double rel = diff / std::max(std::max(an, fn), 1e-30);

    UNSCOPED_INFO(
        "aw=" << aw << " y_sep=" << y_sep << " z_twist=" << z_twist
              << "  ||analytic||=" << an << "  ||fd||=" << fn
              << "  ||diff||=" << diff << "  rel=" << rel);

    CHECK(rel < 1e-4);
}