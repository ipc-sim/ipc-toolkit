#include <ipc/ccd/ccd.hpp>

#ifdef IPC_TOOLKIT_WITH_CORRECT_CCD
#include <tight_inclusion/inclusion_ccd.hpp>
// clang-format off
#warning "TightInclusion does not include point-edge or point-point CCD. Using EVCTCD instead."
// clang-format on
#include <CTCD.h>
#else
#include <CTCD.h>
#endif

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

bool point_point_ccd(
    const Eigen::Vector3d& p0_t0,
    const Eigen::Vector3d& p1_t0,
    const Eigen::Vector3d& p0_t1,
    const Eigen::Vector3d& p1_t1,
    double& toi,
    double conservative_rescaling)
{
    double initial_distance = sqrt(point_point_distance(p0_t0, p1_t0));
    if (CTCD::vertexVertexCTCD(
            p0_t0, p1_t0, p0_t1, p1_t1,
            /*eta=*/(1 - conservative_rescaling) * initial_distance, toi)) {
        if (toi < 1.0e-6) {
            if (CTCD::vertexVertexCTCD(
                    p0_t0, p1_t0, p0_t1, p1_t1, /*eta=*/0, toi)) {
                toi *= conservative_rescaling;
                return true;
            } else {
                return false;
            }
        }
        return true;
    }
    return false;
}

bool point_edge_ccd_2D(
    const Eigen::Vector2d& p_t0,
    const Eigen::Vector2d& e0_t0,
    const Eigen::Vector2d& e1_t0,
    const Eigen::Vector2d& p_t1,
    const Eigen::Vector2d& e0_t1,
    const Eigen::Vector2d& e1_t1,
    double& toi,
    double conservative_rescaling)
{
    Eigen::Vector2d d0 = p_t1 - p_t0, d1 = e0_t1 - e0_t0, d2 = e1_t1 - e1_t0;

    double a = d0[0] * (d2[1] - d1[1]) + d0[1] * (d1[0] - d2[0]) + d2[0] * d1[1]
        - d2[1] * d1[0];
    double b = p_t0[0] * (d2[1] - d1[1]) + d0[0] * (e1_t0[1] - e0_t0[1])
        + d0[1] * (e0_t0[0] - e1_t0[0]) + p_t0[1] * (d1[0] - d2[0])
        + d1[1] * e1_t0[0] + d2[0] * e0_t0[1] - d1[0] * e1_t0[1]
        - d2[1] * e0_t0[0];
    double c = p_t0[0] * (e1_t0[1] - e0_t0[1]) + p_t0[1] * (e0_t0[0] - e1_t0[0])
        + e1_t0[0] * e0_t0[1] - e1_t0[1] * e0_t0[0];

    double roots[2];
    int rootAmt = 0;
    if (a == 0) {
        if (b == 0) {
            // parallel motion, only need to handle colinear case
            if (c == 0) {
                // colinear PP CCD
                if ((p_t0 - e0_t0).dot(d0 - d1) < 0) {
                    roots[rootAmt] = std::sqrt(
                        (p_t0 - e0_t0).squaredNorm() / (d0 - d1).squaredNorm());
                    if (roots[rootAmt] > 0 && roots[rootAmt] <= 1) {
                        ++rootAmt;
                    }
                }
                if ((p_t0 - e1_t0).dot(d0 - d2) < 0) {
                    roots[rootAmt] = std::sqrt(
                        (p_t0 - e1_t0).squaredNorm() / (d0 - d2).squaredNorm());
                    if (roots[rootAmt] > 0 && roots[rootAmt] <= 1) {
                        ++rootAmt;
                    }
                }

                if (rootAmt == 2) {
                    toi = std::min(roots[0], roots[1]) * conservative_rescaling;
                    return true;
                } else if (rootAmt == 1) {
                    toi = roots[0] * conservative_rescaling;
                    return true;
                } else {
                    return false;
                }
            }
        } else {
            rootAmt = 1;
            roots[0] = -c / b;
        }
    } else {
        double delta = b * b - 4 * a * c;
        if (delta == 0) {
            rootAmt = 1;
            roots[0] = -b / (2 * a);
        } else if (delta > 0) {
            rootAmt = 2;
            // accurate expression differs in b's sign
            if (b > 0) {
                roots[0] = (-b - std::sqrt(delta)) / (2 * a);
                roots[1] = 2 * c / (-b - std::sqrt(delta));
            } else {
                roots[0] = 2 * c / (-b + std::sqrt(delta));
                roots[1] = (-b + std::sqrt(delta)) / (2 * a);
            }

            if (roots[0] > roots[1]) {
                std::swap(roots[0], roots[1]);
            }
        }
    }

    for (int i = 0; i < rootAmt; ++i) {
        if (roots[i] > 0 && roots[i] <= 1) {
            // check overlap
            double ratio;
            if (point_edge_distance_type(
                    Eigen::Vector2d(p_t0 + roots[i] * d0),
                    Eigen::Vector2d(e0_t0 + roots[i] * d1),
                    Eigen::Vector2d(e1_t0 + roots[i] * d2))
                == PointEdgeDistanceType::P_E) {
                toi = roots[i] * conservative_rescaling; // TODO: distance eta
                return true;
            }
        }
    }

    return false;
}

bool point_edge_ccd_3D(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& e0_t0,
    const Eigen::Vector3d& e1_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& e0_t1,
    const Eigen::Vector3d& e1_t1,
    double& toi,
    double conservative_rescaling)
{
    double initial_distance = sqrt(point_edge_distance(p_t0, e0_t0, e1_t0));
    if (CTCD::vertexEdgeCTCD(
            p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1,
            /*eta=*/(1 - conservative_rescaling) * initial_distance, toi)) {
        if (toi < 1.0e-6) {
            if (CTCD::vertexEdgeCTCD(
                    p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, /*eta=*/0, toi)) {
                toi *= conservative_rescaling;
                return true;
            } else {
                return false;
            }
        }
        return true;
    }
    return false;
}

bool point_edge_ccd(
    const Eigen::VectorX3d& p_t0,
    const Eigen::VectorX3d& e0_t0,
    const Eigen::VectorX3d& e1_t0,
    const Eigen::VectorX3d& p_t1,
    const Eigen::VectorX3d& e0_t1,
    const Eigen::VectorX3d& e1_t1,
    double& toi,
    double conservative_rescaling)
{
    int dim = p_t0.size();
    assert(e0_t0.size() == dim);
    assert(e1_t0.size() == dim);
    assert(p_t1.size() == dim);
    assert(e0_t1.size() == dim);
    assert(e1_t1.size() == dim);
    if (dim == 2) {
        return point_edge_ccd_2D(
            p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, //
            toi, conservative_rescaling);
    } else {
        return point_edge_ccd_3D(
            p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, //
            toi, conservative_rescaling);
    }
}

bool edge_edge_ccd(
    const Eigen::Vector3d& ea0_t0,
    const Eigen::Vector3d& ea1_t0,
    const Eigen::Vector3d& eb0_t0,
    const Eigen::Vector3d& eb1_t0,
    const Eigen::Vector3d& ea0_t1,
    const Eigen::Vector3d& ea1_t1,
    const Eigen::Vector3d& eb0_t1,
    const Eigen::Vector3d& eb1_t1,
    double& toi,
    double tmax,
    double tolerance,
    int max_iterations,
    double conservative_rescaling)
{
    assert(tmax >= 0 && tmax <= 1.0);

    double initial_distance =
        sqrt(edge_edge_distance(ea0_t0, ea1_t0, eb0_t0, eb1_t0));

    if (initial_distance == 0) {
#ifdef IPC_TOOLKIT_WITH_LOGGER
        logger().warn("Initial edge-edge distance is 0, returning toi=0!");
#endif
        toi = 0;
        return true;
    }

    double min_distance = (1.0 - conservative_rescaling) * initial_distance;

    double output_tolerance = tolerance;
#ifdef IPC_TOOLKIT_WITH_CORRECT_CCD
    // 0: normal ccd method which only checks t = [0,1]
    // 1: ccd with max_itr and t=[0, t_max]
    const int CCD_TYPE = 1;
    bool is_impacting = inclusion_ccd::edgeEdgeCCD_double(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
        { { -1, -1, -1 } }, // rounding error (auto)
        min_distance,       // minimum separation distance
        toi,                // time of impact
        tolerance,          // delta
        tmax,               // maximum time to check
        max_iterations,     // maximum number of iterations
        output_tolerance,   // delta_actual
        CCD_TYPE);
#else
    bool is_impacting = CTCD::edgeEdgeCTCD(
        ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
        min_distance, toi);
#endif

    if (is_impacting && toi < 1.0e-6) {
#ifdef IPC_TOOLKIT_WITH_CORRECT_CCD
        is_impacting = inclusion_ccd::edgeEdgeCCD_double(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
            { { -1, -1, -1 } }, // rounding error (auto)
            0,                  // minimum separation distance
            toi,                // time of impact
            tolerance,          // delta
            tmax,               // maximum time to check
            max_iterations,     // maximum number of iterations
            output_tolerance,   // delta_actual
            CCD_TYPE);
#else
        is_impacting = CTCD::edgeEdgeCTCD(
            ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
            /*eta=*/0, toi);
#endif

        if (is_impacting) {
            toi *= conservative_rescaling;
#ifdef IPC_TOOLKIT_WITH_LOGGER
            if (toi == 0) {
                logger().warn(
                    "Edge-edge CCD is overly conservative "
                    "(initial_distance={:g}, actual_tolerance={:g})!",
                    initial_distance, output_tolerance);
            }
#endif
        }
    }

    return is_impacting;
}

bool point_triangle_ccd(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& t0_t0,
    const Eigen::Vector3d& t1_t0,
    const Eigen::Vector3d& t2_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& t0_t1,
    const Eigen::Vector3d& t1_t1,
    const Eigen::Vector3d& t2_t1,
    double& toi,
    double tmax,
    double tolerance,
    int max_iterations,
    double conservative_rescaling)
{
    assert(tmax >= 0 && tmax <= 1.0);

    double initial_distance =
        sqrt(point_triangle_distance(p_t0, t0_t0, t1_t0, t2_t0));

    if (initial_distance == 0) {
#ifdef IPC_TOOLKIT_WITH_LOGGER
        logger().warn("Initial point-triangle distance is 0, returning toi=0!");
#endif
        toi = 0;
        return true;
    }

    double min_distance = (1.0 - conservative_rescaling) * initial_distance;

    double output_tolerance = tolerance;
#ifdef IPC_TOOLKIT_WITH_CORRECT_CCD
    // 0: normal ccd method which only checks t = [0,1]
    // 1: ccd with max_itr and t=[0, t_max]
    const int CCD_TYPE = 1;
    bool is_impacting = inclusion_ccd::vertexFaceCCD_double(
        p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1,
        { { -1, -1, -1 } }, // rounding error (auto)
        min_distance,       // minimum separation distance
        toi,                // time of impact
        tolerance,          // delta
        tmax,               // maximum time to check
        max_iterations,     // maximum number of iterations
        output_tolerance,   // delta_actual
        CCD_TYPE);
#else
    bool is_impacting = CTCD::vertexFaceCTCD(
        p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1, //
        min_distance, toi);
#endif

    if (is_impacting && toi < 1e-6) {

#ifdef IPC_TOOLKIT_WITH_CORRECT_CCD
        is_impacting = inclusion_ccd::vertexFaceCCD_double(
            p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1,
            { { -1, -1, -1 } }, // rounding error (auto)
            0,                  // minimum separation distance
            toi,                // time of impact
            tolerance,          // delta
            tmax,               // maximum time to check
            max_iterations,     // maximum number of iterations
            output_tolerance,   // delta_actual
            CCD_TYPE);
#else
        is_impacting = CTCD::vertexFaceCTCD(
            p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1,
            /*eta=*/0, toi);
#endif

        if (is_impacting) {
            toi *= conservative_rescaling;
#ifdef IPC_TOOLKIT_WITH_LOGGER
            if (toi == 0) {
                logger().warn(
                    "Point-triangle CCD is overly conservative "
                    "(initial_distance={:g}, actual_tolerance={:g})!",
                    initial_distance, output_tolerance);
            }
#endif
        }
    }

    return is_impacting;
}

} // namespace ipc
