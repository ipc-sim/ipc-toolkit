#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <tests/utils.hpp>

#include <ipc/distance/point_edge.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <finitediff.hpp>

using namespace ipc;

namespace {
template <int dim> double point_edge_distance_stacked(const Eigen::VectorXd& x)
{
    assert(x.size() == 3 * dim);
    return point_edge_distance(
        x.head<dim>(), x.segment<dim>(dim), x.tail<dim>());
}
} // namespace

TEMPLATE_TEST_CASE_SIG(
    "Point-edge distance",
    "[distance][point-edge][gradient][hessian]",
    ((int dim), dim),
    2,
    3)
{
    double expected_distance = GENERATE(-10, -1, -1e-12, 0, 1e-12, 1, 10);
    VectorMax3d p = VectorMax3d::Zero(dim);
    p.y() = expected_distance;
    VectorMax3d e0 = VectorMax3d::Zero(dim);
    e0.x() = -10;
    VectorMax3d e1 = VectorMax3d::Zero(dim);
    e1.x() = 10;

    { // Distance
        double distance = point_edge_distance(p, e0, e1);
        CHECK(distance == Catch::Approx(expected_distance * expected_distance));
    }

    { //  Gradient
        const VectorMax9d grad = point_edge_distance_gradient(p, e0, e1);

        // Compute the gradient using finite differences
        VectorMax9d x(3 * dim);
        x << p, e0, e1;
        Eigen::VectorXd fgrad;
        fd::finite_gradient(x, point_edge_distance_stacked<dim>, fgrad);

        CHECK(fd::compare_gradient(grad, fgrad));
    }

    { // Hessian
        const MatrixMax9d hess = point_edge_distance_hessian(p, e0, e1);

        // Compute the gradient using finite differences
        VectorMax9d x(3 * dim);
        x << p, e0, e1;
        Eigen::MatrixXd fhess;
        fd::finite_hessian(x, point_edge_distance_stacked<dim>, fhess);

        CHECK(fd::compare_hessian(hess, fhess, 1e-2));
    }
}

TEMPLATE_TEST_CASE_SIG(
    "Point-edge distance all types (random edges)",
    "[distance][point-edge][gradient][hessian]",
    ((int dim), dim),
    2,
    3)
{
    const double alpha = GENERATE(range(-1.0, 2.0, 0.1));
    const double d = GENERATE(range(-10.0, 10.0, 1.0));
    const int n_random_edges = 20;

    for (int i = 0; i < n_random_edges; i++) {
        VectorMax3d e0, e1, n;
        if (dim == 2) {
            e0 = Eigen::Vector2d::Random();
            e1 = Eigen::Vector2d::Random();
            n = tests::edge_normal(e0, e1);
        } else {
            e0 = Eigen::Vector3d::Random();
            e1 = Eigen::Vector3d::Random();
            n = Eigen::Vector3d(e1 - e0).cross(Eigen::Vector3d::UnitX());
            n.normalize();
        }

        const VectorMax3d p = ((e1 - e0) * alpha + e0) + d * n;

        CAPTURE(alpha, dim);

        { // Distance
            double expected_distance =
                alpha < 0 ? (e0 - p).norm() : (alpha > 1 ? (e1 - p).norm() : d);
            expected_distance *= expected_distance;
            const double distance = point_edge_distance(p, e0, e1);
            CHECK(distance == Catch::Approx(expected_distance).margin(1e-15));
        }

        // Gradient (skip C1 transition points)
        if (abs(alpha) < 1e-5 && abs(alpha - 1.0) < 1e-5) {
            const VectorMax9d grad = point_edge_distance_gradient(p, e0, e1);

            // Compute the gradient using finite differences
            VectorMax9d x(3 * dim);
            x << p, e0, e1;

            Eigen::VectorXd fgrad;
            fd::finite_gradient(x, point_edge_distance_stacked<dim>, fgrad);

            CHECK(fd::compare_gradient(grad, fgrad));
        }

        // Hessian (skip C1 transition points)
        if (abs(alpha) < 1e-5 && abs(alpha - 1.0) < 1e-5) {
            const MatrixMax9d hess = point_edge_distance_hessian(p, e0, e1);
            // Compute the gradient using finite differences
            VectorMax9d x(3 * dim);
            x << p, e0, e1;
            Eigen::MatrixXd fhess;
            fd::finite_hessian(x, point_edge_distance_stacked<dim>, fhess);
            CHECK(fd::compare_hessian(hess, fhess, 1e-2));
        }
    }
}

TEMPLATE_TEST_CASE_SIG(
    "Point-edge distance all types",
    "[distance][point-edge][gradient][hessian]",
    ((int dim), dim),
    2,
    3)
{
    const double alpha = GENERATE(range(-1.0, 2.0, 0.1));
    const double d = GENERATE(range(-10.0, 10.0, 1.0));

    VectorMax3d e0 = VectorMax3d::Zero(dim);
    e0.x() = -10;
    VectorMax3d e1 = VectorMax3d::Zero(dim);
    e1.x() = 10;

    const VectorMax3d p =
        ((e1 - e0) * alpha + e0) + d * Eigen::Vector3d::UnitY().head(dim);

    CAPTURE(alpha, d);

    { // Distance type
        const auto dtype = point_edge_distance_type(p, e0, e1);
        if (alpha < -1e-8) {
            CHECK(dtype == PointEdgeDistanceType::P_E0);
        } else if (alpha > 1 + 1e-8) {
            CHECK(dtype == PointEdgeDistanceType::P_E1);
        } else if (alpha > 1e-8 && alpha < 1 - 1e-8) {
            CHECK(dtype == PointEdgeDistanceType::P_E);
        }
        // ignore alpha ≈ 0 and alpha ≈ 1
    }

    { // Distance
        const double expected_distance = alpha < 0
            ? (e0 - p).squaredNorm()
            : (alpha > 1 ? (e1 - p).squaredNorm() : (d * d));
        const double distance = point_edge_distance(p, e0, e1);
        CHECK(distance == Catch::Approx(expected_distance).margin(1e-15));
    }

    // Gradient (skip C1 transition points)
    // if (abs(alpha) < 1e-5 && abs(alpha - 1.0) < 1e-5) {
    {
        const VectorMax9d grad = point_edge_distance_gradient(p, e0, e1);

        // Compute the gradient using finite differences
        VectorMax9d x(3 * dim);
        x << p, e0, e1;
        Eigen::VectorXd fgrad;
        fd::finite_gradient(x, point_edge_distance_stacked<dim>, fgrad);

        CHECK(fd::compare_gradient(grad, fgrad));
    }

    // Hessian (skip C1 transition points)
    if (abs(alpha) < 1e-5 && abs(alpha - 1.0) < 1e-5) {
        const MatrixMax9d hess = point_edge_distance_hessian(p, e0, e1);
        // Compute the gradient using finite differences
        VectorMax9d x(3 * dim);
        x << p, e0, e1;
        Eigen::MatrixXd fhess;
        fd::finite_hessian(x, point_edge_distance_stacked<dim>, fhess);
        CHECK(fd::compare_hessian(hess, fhess, 1e-2));
    }
}
