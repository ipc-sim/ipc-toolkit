#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>
#include <catch2/generators/catch_generators_random.hpp>

#include <ipc/distance/point_plane.hpp>

#include <finitediff.hpp>

using namespace ipc;

namespace {
double point_plane_distance_stacked(const Eigen::VectorXd& x)
{
    assert(x.size() == 12);
    return point_plane_distance(
        x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>());
};
} // namespace

TEST_CASE(
    "Point-plane distance and derivatives (dynamic plane)",
    "[distance][point-plane][gradient][hessian]")
{
    double x = GENERATE(take(10, random(-10.0, 10.0)));
    double y = GENERATE(take(10, random(-10.0, 10.0)));
    double z = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d p(x, y, z);

    double y_plane = GENERATE(take(10, random(-10.0, 10.0)));
    Eigen::Vector3d t0(-1, y_plane, 0), t1(1, y_plane, -1), t2(1, y_plane, 0);

    {
        double distance = point_plane_distance(p, t0, t1, t2);
        double expected_distance = std::abs(y - y_plane);
        CHECK(distance == Catch::Approx(expected_distance * expected_distance));
    }

    {
        const Vector12d grad = point_plane_distance_gradient(p, t0, t1, t2);

        Vector12d x_vec;
        x_vec << p, t0, t1, t2;
        Eigen::VectorXd expected_grad;
        expected_grad.resize(grad.size());
        fd::finite_gradient(x_vec, point_plane_distance_stacked, expected_grad);

        CAPTURE((grad - expected_grad).norm());
        CHECK(fd::compare_gradient(grad, expected_grad));
    }

    {
        const Matrix12d hess = point_plane_distance_hessian(p, t0, t1, t2);

        Vector12d x_vec;
        x_vec << p, t0, t1, t2;
        Eigen::MatrixXd expected_hess;
        expected_hess.resize(hess.rows(), hess.cols());
        fd::finite_hessian(x_vec, point_plane_distance_stacked, expected_hess);

        CAPTURE((hess - expected_hess).norm());
        CHECK(fd::compare_hessian(hess, expected_hess, 5e-2));
    }
}

TEST_CASE(
    "Point-plane distance and derivatives (static plane)",
    "[distance][point-plane][implicit]")
{
    const Eigen::Vector3d p(0, 2, 0);
    const Eigen::Vector3d origin(0, 0, 0);
    const Eigen::Vector3d normal(0, 1, 0);

    CHECK(point_plane_distance(p, origin, normal) == p.y() * p.y());
    CHECK(point_plane_distance_gradient(p, origin, normal) == 2 * p);
    {
        Eigen::Matrix3d Ha = point_plane_distance_hessian(p, origin, normal);
        Eigen::Matrix3d He = Eigen::Matrix3d::Zero();
        He(1, 1) = 2;
        CHECK(Ha == He);
    }
}