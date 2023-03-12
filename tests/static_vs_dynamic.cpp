#include <catch2/catch_all.hpp>

#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/eigen_ext.hpp>

using namespace ipc;

TEST_CASE("Template dynamic vs static", "[!benchmark][eigen]")
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Random(100, 3);

    int vi = 0, t0i = 50, t1i = 75, t2i = 99;

    BENCHMARK("Static")
    {
        MatrixMax12d hess;

        Eigen::Vector3d p = V.row(vi);
        Eigen::Vector3d t0 = V.row(t0i);
        Eigen::Vector3d t1 = V.row(t1i);
        Eigen::Vector3d t2 = V.row(t2i);

        point_triangle_distance_hessian(p, t0, t1, t2, hess);
    };

    BENCHMARK("Dynamic")
    {
        MatrixMax12d hess;
        point_triangle_distance_hessian(
            V.row(vi), V.row(t0i), V.row(t1i), V.row(t2i), hess);
    };
}
