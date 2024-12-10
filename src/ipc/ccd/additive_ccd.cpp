//
// Source modified from https://github.com/ipc-sim/Codim-IPC
// under Appache-2.0 License.
//
// Modifications:
//  • remove broad phase functions
//  • refactor code to use a single implementation of the additive_ccd algorithm
//  • utilize our distance function rather than porting the Codim-IPC versions
//  • return true if the initial distance is less than the minimum distance
//  • add an explicit tmax parameter rather than relying on the initial value of
//    toi
//  • add a maximum number of iterations to limit the computation time
//
// NOTE: These methods are provided for reference comparison with [Li et al.
// 2021] and is not utilized by the high-level functionality. In compairson to
// Tight Inclusion CCD, this CCD method is not provably conservative and so can
// potentially produce false negatives (i.e., miss collisions) due to
// floating-point rounding error. However, it is much faster than Tight
// Inclusion CCD (>100×) and very robust due to the gaps and conservative
// rescaling used.
//

#include "additive_ccd.hpp"

#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

namespace {
    template <typename... Args> void subtract_mean(Args&... args)
    {
        constexpr double n = sizeof...(args);
        static_assert(n > 0, "At least one argument is required");

        using T = typename std::tuple_element<0, std::tuple<Args...>>::type;
        const int dim = std::get<0>(std::tuple<const Args&...>(args...)).size();

        T mean = T::Zero(dim);
        for (const T& value : { args... }) {
            mean += value;
        }
        mean /= n;

        for (T* value : { &args... }) {
            (*value) -= mean;
        }
    }

    VectorMax12d stack(const VectorMax3d& x) { return x; }

    template <typename... Args>
    VectorMax12d stack(const VectorMax3d& x0, const Args&... args)
    {
        VectorMax12d x(x0.size() * (1 + sizeof...(args)));
        x.head(x0.size()) = x0;
        x.tail(x0.size() * sizeof...(args)) = stack(args...);
        return x;
    }
} // namespace

AdditiveCCD::AdditiveCCD(
    const long _max_iterations, const double _conservative_rescaling)
    : max_iterations(_max_iterations)
    , conservative_rescaling(_conservative_rescaling)
{
}

bool AdditiveCCD::additive_ccd(
    VectorMax12d x,
    const VectorMax12d& dx,
    const std::function<double(const VectorMax12d&)>& distance_squared,
    const double max_disp_mag,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    assert(conservative_rescaling > 0 && conservative_rescaling <= 1);

    const double min_distance_sq = min_distance * min_distance;

    double d, d_sq;
    d = std::sqrt(d_sq = distance_squared(x));
    assert(d > min_distance);

    double d_func = d_sq - min_distance_sq;
    assert(d_func > 0);
    const double gap = // (d - ξ) = (d² - ξ²) / (d + ξ)
        (1 - conservative_rescaling) * d_func / (d + min_distance);
    if (gap < std::numeric_limits<double>::epsilon()) {
        logger().warn(
            "Small gap {:g} ≤ ε in Additive CCD can lead to missed collisions",
            gap);
    }

    toi = 0;
    for (long i = 0; max_iterations < 0 || i < max_iterations; ++i) {
        // tₗ = η ⋅ (d - ξ) / lₚ = η ⋅ (d² - ξ²) / (lₚ ⋅ (d + ξ))
        const double toi_lower_bound = conservative_rescaling * d_func
            / ((d + min_distance) * max_disp_mag);

        x += toi_lower_bound * dx;

        d = std::sqrt(d_sq = distance_squared(x));

        d_func = d_sq - min_distance_sq;
        assert(d_func > 0);
        if (toi > 0 && d_func / (d + min_distance) < gap) {
            break; // distance (including thickness) is less than gap
        }

        toi += toi_lower_bound;
        if (toi > tmax) {
            return false; // collision occurs after tmax
        }

        if (max_iterations < 0 && i == DEFAULT_MAX_ITERATIONS) {
            logger().warn(
                "Slow convergence in Additive CCD. Perhaps the gap is too small (gap={:g})?",
                gap);
        }
    }

    return true;
}

bool AdditiveCCD::point_point_ccd(
    const VectorMax3d& p0_t0,
    const VectorMax3d& p1_t0,
    const VectorMax3d& p0_t1,
    const VectorMax3d& p1_t1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    const int dim = p0_t0.size();
    assert(dim == p1_t0.size() && dim == p0_t1.size() && dim == p1_t1.size());

    const double initial_distance = point_point_distance(p0_t0, p1_t0);
    if (initial_distance <= min_distance * min_distance) {
        logger().warn(
            "Initial distance {} ≤ d_min={}, returning toi=0!",
            std::sqrt(initial_distance), min_distance);
        toi = 0;
        return true;
    }

    VectorMax3d dp0 = p0_t1 - p0_t0;
    VectorMax3d dp1 = p1_t1 - p1_t0;
    subtract_mean(dp0, dp1);

    const double max_disp_mag = dp0.norm() + dp1.norm();
    if (max_disp_mag == 0) {
        return false;
    }

    auto distance_squared = [dim](const VectorMax12d& x) {
        return point_point_distance(x.head(dim), x.tail(dim));
    };

    const VectorMax12d x = stack(p0_t0, p1_t0);
    const VectorMax12d dx = stack(dp0, dp1);

    return additive_ccd(
        x, dx, distance_squared, max_disp_mag, toi, min_distance, tmax);
}

bool AdditiveCCD::point_edge_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& e0_t0,
    const VectorMax3d& e1_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& e0_t1,
    const VectorMax3d& e1_t1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    const int dim = p_t0.size();
    assert(dim == e0_t0.size() && dim == e1_t0.size());
    assert(dim == p_t1.size() && dim == e0_t1.size() && dim == e1_t1.size());

    const double initial_distance = point_edge_distance(p_t0, e0_t0, e1_t0);
    if (initial_distance <= min_distance * min_distance) {
        logger().warn(
            "Initial distance {} ≤ d_min={}, returning toi=0!",
            std::sqrt(initial_distance), min_distance);
        toi = 0;
        return true;
    }

    VectorMax3d dp = p_t1 - p_t0;
    VectorMax3d de0 = e0_t1 - e0_t0;
    VectorMax3d de1 = e1_t1 - e1_t0;
    subtract_mean(dp, de0, de1);

    const double max_disp_mag =
        dp.norm() + std::sqrt(std::max(de0.squaredNorm(), de1.squaredNorm()));
    if (max_disp_mag == 0) {
        return false;
    }

    auto distance_squared = [dim](const VectorMax12d& x) {
        return point_edge_distance(
            x.head(dim), x.segment(dim, dim), x.tail(dim));
    };

    const VectorMax12d x = stack(p_t0, e0_t0, e1_t0);
    const VectorMax12d dx = stack(dp, de0, de1);

    return additive_ccd(
        x, dx, distance_squared, max_disp_mag, toi, min_distance, tmax);
}

bool AdditiveCCD::point_triangle_ccd(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& t0_t0,
    const Eigen::Vector3d& t1_t0,
    const Eigen::Vector3d& t2_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& t0_t1,
    const Eigen::Vector3d& t1_t1,
    const Eigen::Vector3d& t2_t1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    const double initial_distance =
        point_triangle_distance(p_t0, t0_t0, t1_t0, t2_t0);
    if (initial_distance <= min_distance * min_distance) {
        logger().warn(
            "Initial distance {} ≤ d_min={}, returning toi=0!",
            std::sqrt(initial_distance), min_distance);
        toi = 0;
        return true;
    }

    Eigen::Vector3d dp = p_t1 - p_t0;
    Eigen::Vector3d dt0 = t0_t1 - t0_t0;
    Eigen::Vector3d dt1 = t1_t1 - t1_t0;
    Eigen::Vector3d dt2 = t2_t1 - t2_t0;
    subtract_mean(dp, dt0, dt1, dt2);

    const double max_disp_mag = dp.norm()
        + std::sqrt(std::max(
            { dt0.squaredNorm(), dt1.squaredNorm(), dt2.squaredNorm() }));
    if (max_disp_mag == 0) {
        return false;
    }

    auto distance_squared = [](const VectorMax12d& x) {
        return point_triangle_distance(
            x.head<3>(), x.segment<3>(3), x.segment<3>(6), x.tail<3>());
    };

    const VectorMax12d x = stack(p_t0, t0_t0, t1_t0, t2_t0);
    const VectorMax12d dx = stack(dp, dt0, dt1, dt2);

    return additive_ccd(
        x, dx, distance_squared, max_disp_mag, toi, min_distance, tmax);
}

bool AdditiveCCD::edge_edge_ccd(
    const Eigen::Vector3d& ea0_t0,
    const Eigen::Vector3d& ea1_t0,
    const Eigen::Vector3d& eb0_t0,
    const Eigen::Vector3d& eb1_t0,
    const Eigen::Vector3d& ea0_t1,
    const Eigen::Vector3d& ea1_t1,
    const Eigen::Vector3d& eb0_t1,
    const Eigen::Vector3d& eb1_t1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    const double initial_distance =
        edge_edge_distance(ea0_t0, ea1_t0, eb0_t0, eb1_t0);
    if (initial_distance <= min_distance * min_distance) {
        logger().warn(
            "Initial distance {} ≤ d_min={}, returning toi=0!",
            std::sqrt(initial_distance), min_distance);
        toi = 0;
        return true;
    }

    Eigen::Vector3d dea0 = ea0_t1 - ea0_t0;
    Eigen::Vector3d dea1 = ea1_t1 - ea1_t0;
    Eigen::Vector3d deb0 = eb0_t1 - eb0_t0;
    Eigen::Vector3d deb1 = eb1_t1 - eb1_t0;
    subtract_mean(dea0, dea1, deb0, deb1);

    const double max_disp_mag =
        std::sqrt(std::max(dea0.squaredNorm(), dea1.squaredNorm()))
        + std::sqrt(std::max(deb0.squaredNorm(), deb1.squaredNorm()));
    if (max_disp_mag == 0) {
        return false;
    }

    const double min_distance_sq = min_distance * min_distance;
    auto distance_squared = [min_distance_sq](const VectorMax12d& x) {
        const auto& ea0 = x.head<3>();
        const auto& ea1 = x.segment<3>(3);
        const auto& eb0 = x.segment<3>(6);
        const auto& eb1 = x.tail<3>();

        double d_sq = edge_edge_distance(ea0, ea1, eb0, eb1);
        if (d_sq - min_distance_sq <= 0) {
            // since we ensured other place that all dist smaller than d̂ are
            // positive, this must be some far away nearly parallel edges
            d_sq = std::min(
                { (ea0 - eb0).squaredNorm(), (ea0 - eb1).squaredNorm(),
                  (ea1 - eb0).squaredNorm(), (ea1 - eb1).squaredNorm() });
        }
        return d_sq;
    };

    const VectorMax12d x = stack(ea0_t0, ea1_t0, eb0_t0, eb1_t0);
    const VectorMax12d dx = stack(dea0, dea1, deb0, deb1);

    return additive_ccd(
        x, dx, distance_squared, max_disp_mag, toi, min_distance, tmax);
}

} // namespace ipc