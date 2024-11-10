#include "intersection.hpp"

#include <ipc/config.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Geometry>
#include <igl/predicates/predicates.h>

#ifdef IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION
#include <rational/rational.hpp>
#endif

namespace ipc {

#ifdef IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION
namespace {
    bool is_edge_intersecting_triangle_rational(
        const Eigen::Vector3d& e0_float,
        const Eigen::Vector3d& e1_float,
        const Eigen::Vector3d& t0_float,
        const Eigen::Vector3d& t1_float,
        const Eigen::Vector3d& t2_float)
    {
        using namespace rational;

        typedef Eigen::Matrix<
            Rational, 3, 1, Eigen::ColMajor | Eigen::DontAlign>
            Vector3r;

        Vector3r e0, e1, t0, t1, t2;

        for (int d = 0; d < 3; ++d) {
            e0[d] = e0_float[d];
            e1[d] = e1_float[d];

            t0[d] = t0_float[d];
            t1[d] = t1_float[d];
            t2[d] = t2_float[d];
        }

        const Rational d = e0[0] * t0[1] * t1[2] - e0[0] * t0[1] * t2[2]
            - e0[0] * t0[2] * t1[1] + e0[0] * t0[2] * t2[1]
            + e0[0] * t1[1] * t2[2] - e0[0] * t1[2] * t2[1]
            - e0[1] * t0[0] * t1[2] + e0[1] * t0[0] * t2[2]
            + e0[1] * t0[2] * t1[0] - e0[1] * t0[2] * t2[0]
            - e0[1] * t1[0] * t2[2] + e0[1] * t1[2] * t2[0]
            + e0[2] * t0[0] * t1[1] - e0[2] * t0[0] * t2[1]
            - e0[2] * t0[1] * t1[0] + e0[2] * t0[1] * t2[0]
            + e0[2] * t1[0] * t2[1] - e0[2] * t1[1] * t2[0]
            - e1[0] * t0[1] * t1[2] + e1[0] * t0[1] * t2[2]
            + e1[0] * t0[2] * t1[1] - e1[0] * t0[2] * t2[1]
            - e1[0] * t1[1] * t2[2] + e1[0] * t1[2] * t2[1]
            + e1[1] * t0[0] * t1[2] - e1[1] * t0[0] * t2[2]
            - e1[1] * t0[2] * t1[0] + e1[1] * t0[2] * t2[0]
            + e1[1] * t1[0] * t2[2] - e1[1] * t1[2] * t2[0]
            - e1[2] * t0[0] * t1[1] + e1[2] * t0[0] * t2[1]
            + e1[2] * t0[1] * t1[0] - e1[2] * t0[1] * t2[0]
            - e1[2] * t1[0] * t2[1] + e1[2] * t1[1] * t2[0];
        if (d.sign() == 0) {
            return true;
        }

        // t is the parametric coordinate for the edge
        const Rational t = (e0[0] * t0[1] * t1[2] - e0[0] * t0[1] * t2[2]
                            - e0[0] * t0[2] * t1[1] + e0[0] * t0[2] * t2[1]
                            + e0[0] * t1[1] * t2[2] - e0[0] * t1[2] * t2[1]
                            - e0[1] * t0[0] * t1[2] + e0[1] * t0[0] * t2[2]
                            + e0[1] * t0[2] * t1[0] - e0[1] * t0[2] * t2[0]
                            - e0[1] * t1[0] * t2[2] + e0[1] * t1[2] * t2[0]
                            + e0[2] * t0[0] * t1[1] - e0[2] * t0[0] * t2[1]
                            - e0[2] * t0[1] * t1[0] + e0[2] * t0[1] * t2[0]
                            + e0[2] * t1[0] * t2[1] - e0[2] * t1[1] * t2[0]
                            - t0[0] * t1[1] * t2[2] + t0[0] * t1[2] * t2[1]
                            + t0[1] * t1[0] * t2[2] - t0[1] * t1[2] * t2[0]
                            - t0[2] * t1[0] * t2[1] + t0[2] * t1[1] * t2[0])
            / d;

        if (t < 0 || t > 1) {
            return false;
        }

        // u is the first barycentric coordinate for the triangle
        const Rational u = (-e0[0] * e1[1] * t0[2] + e0[0] * e1[1] * t2[2]
                            + e0[0] * e1[2] * t0[1] - e0[0] * e1[2] * t2[1]
                            - e0[0] * t0[1] * t2[2] + e0[0] * t0[2] * t2[1]
                            + e0[1] * e1[0] * t0[2] - e0[1] * e1[0] * t2[2]
                            - e0[1] * e1[2] * t0[0] + e0[1] * e1[2] * t2[0]
                            + e0[1] * t0[0] * t2[2] - e0[1] * t0[2] * t2[0]
                            - e0[2] * e1[0] * t0[1] + e0[2] * e1[0] * t2[1]
                            + e0[2] * e1[1] * t0[0] - e0[2] * e1[1] * t2[0]
                            - e0[2] * t0[0] * t2[1] + e0[2] * t0[1] * t2[0]
                            + e1[0] * t0[1] * t2[2] - e1[0] * t0[2] * t2[1]
                            - e1[1] * t0[0] * t2[2] + e1[1] * t0[2] * t2[0]
                            + e1[2] * t0[0] * t2[1] - e1[2] * t0[1] * t2[0])
            / d;
        // v is the second barycentric coordinate for the triangle
        const Rational v = (e0[0] * e1[1] * t0[2] - e0[0] * e1[1] * t1[2]
                            - e0[0] * e1[2] * t0[1] + e0[0] * e1[2] * t1[1]
                            + e0[0] * t0[1] * t1[2] - e0[0] * t0[2] * t1[1]
                            - e0[1] * e1[0] * t0[2] + e0[1] * e1[0] * t1[2]
                            + e0[1] * e1[2] * t0[0] - e0[1] * e1[2] * t1[0]
                            - e0[1] * t0[0] * t1[2] + e0[1] * t0[2] * t1[0]
                            + e0[2] * e1[0] * t0[1] - e0[2] * e1[0] * t1[1]
                            - e0[2] * e1[1] * t0[0] + e0[2] * e1[1] * t1[0]
                            + e0[2] * t0[0] * t1[1] - e0[2] * t0[1] * t1[0]
                            - e1[0] * t0[1] * t1[2] + e1[0] * t0[2] * t1[1]
                            + e1[1] * t0[0] * t1[2] - e1[1] * t0[2] * t1[0]
                            - e1[2] * t0[0] * t1[1] + e1[2] * t0[1] * t1[0])
            / d;

        return u >= 0 && u <= 1 && v >= 0 && v <= 1 && u + v <= 1;
    }
} // namespace
#endif

bool is_edge_intersecting_triangle(
    const Eigen::Vector3d& e0,
    const Eigen::Vector3d& e1,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2)
{
    igl::predicates::exactinit();
    const auto ori1 = igl::predicates::orient3d(t0, t1, t2, e0);
    const auto ori2 = igl::predicates::orient3d(t0, t1, t2, e1);

    if (ori1 != igl::predicates::Orientation::COPLANAR
        && ori2 != igl::predicates::Orientation::COPLANAR && ori1 == ori2) {
        // edge is completly on one side of the plane that triangle is in
        return false;
    }

#ifdef IPC_TOOLKIT_WITH_RATIONAL_INTERSECTION
    return is_edge_intersecting_triangle_rational(e0, e1, t0, t1, t2);
#else
    Eigen::Matrix3d M;
    M.col(0) = t1 - t0;
    M.col(1) = t2 - t0;
    M.col(2) = e0 - e1;
    Eigen::Vector3d uvt = M.fullPivLu().solve(e0 - t0);
    return uvt[0] >= 0.0 && uvt[1] >= 0.0 && uvt[0] + uvt[1] <= 1.0
        && uvt[2] >= 0.0 && uvt[2] <= 1.0;
#endif
}

} // namespace ipc
