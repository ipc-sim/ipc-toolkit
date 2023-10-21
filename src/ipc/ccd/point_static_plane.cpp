#include "point_static_plane.hpp"

#include <ipc/distance/point_plane.hpp>

#include <array>

namespace ipc {

inline bool is_in_01(double x) { return 0 <= x && x <= 1; };

bool point_static_plane_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& plane_origin,
    const VectorMax3d& plane_normal,
    double& toi,
    double conservative_rescaling)
{
    static constexpr double SMALL_TOI = 1e-6;

    assert(p_t1.size() == p_t0.size());
    assert(plane_origin.size() == p_t0.size());
    assert(plane_normal.size() == p_t0.size());

    double initial_distance =
        sqrt(point_plane_distance(p_t0, plane_origin, plane_normal));

    if (initial_distance == 0) {
        logger().warn("Initial point-plane distance is 0, returning toi=0!");
        toi = 0;
        return true;
    }

    auto compute_toi = [&](double d) -> double {
        return (d * plane_normal.norm() + plane_normal.dot(plane_origin - p_t0))
            / plane_normal.dot(p_t1 - p_t0);
    };
    auto compute_tois = [&compute_toi](double d) -> std::array<double, 2> {
        return { { compute_toi(d), compute_toi(-d) } };
    };

    double min_distance = (1.0 - conservative_rescaling) * initial_distance;
    assert(min_distance < initial_distance);
    std::array<double, 2> tois = compute_tois(min_distance);

    bool is_impacting = is_in_01(tois[0]) || is_in_01(tois[1]);
    if (is_in_01(tois[0]) && is_in_01(tois[1])) {
        toi = std::min(tois[0], tois[1]);
    } else if (is_in_01(tois[0])) {
        toi = tois[0];
    } else if (is_in_01(tois[1])) {
        toi = tois[1];
    }

    if (is_impacting && toi < SMALL_TOI) {
        toi = compute_toi(/*d=*/0);
        is_impacting = is_in_01(toi);
        if (is_impacting) {
            toi *= conservative_rescaling;
            if (toi == 0) {
                logger().warn(
                    "Point-static plane CCD is overly conservative (toi={:g} "
                    "and trajectory_length={:g}, but initial_distance={:g})!",
                    toi, (p_t1 - p_t0).norm(), initial_distance);
            }
        }
    }

    return is_impacting;
}

} // namespace ipc
