#include "nonlinear_ccd.hpp"

#include <ipc/ccd/check_initial_distance.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <stack>

// #define USE_FIXED_PIECES

namespace ipc {

#ifdef USE_FIXED_PIECES
static constexpr size_t FIXED_NUM_PIECES = 100L;
#else
static constexpr size_t MAX_NUM_SUBDIVISIONS = 1000L;
#endif

// ============================================================================

#ifdef IPC_TOOLKIT_WITH_FILIB
double IntervalNonlinearTrajectory::max_distance_from_linear(
    const double t0, const double t1) const
{
    // Estimate max t ∈ [0, 1] ‖ p((t1 - t0) * t + t0) - lerp(p(t0), p(t1), t) ‖
    const VectorMax3I p_t0 = (*this)(t0).cast<filib::Interval>();
    const VectorMax3I p_t1 = (*this)(t1).cast<filib::Interval>();

    double max_d = 0;

    constexpr double n = 100.0;
    double ti0 = t0;
    for (int i = 1; i <= n; i++) {
        const double ti1 = i / n * (t1 - t0) + t0;

        const VectorMax3I p = (*this)(filib::Interval(ti0, ti1));
        const filib::Interval d = norm(
            p - ((p_t1 - p_t0) * filib::Interval(i / n, (i + 1) / n) + p_t0));

        max_d = std::max(max_d, d.SUP);
        ti0 = ti1;
    }

    return max_d;
}
#endif

// ============================================================================

NonlinearCCD::NonlinearCCD(
    const double _conservative_rescaling,
    std::shared_ptr<NarrowPhaseCCD> _narrow_phase_ccd)
    : conservative_rescaling(_conservative_rescaling)
{
    set_narrow_phase_ccd(std::move(_narrow_phase_ccd));
}

void NonlinearCCD::set_narrow_phase_ccd(
    std::shared_ptr<NarrowPhaseCCD> _narrow_phase_ccd)
{
    if (_narrow_phase_ccd == nullptr) {
        _narrow_phase_ccd = std::make_shared<TightInclusionCCD>();
    }
    // Tell the narrow phase it is a component of this piecewise driver: it
    // must skip its standalone conservative strategy and return the raw
    // per-piece time of impact (see NarrowPhaseCCD::component_mode).
    _narrow_phase_ccd->component_mode = true;
    m_narrow_phase_ccd = std::move(_narrow_phase_ccd);
}

bool NonlinearCCD::conservative_piecewise_linear_ccd(
    const std::function<double(const double)>& distance,
    const std::function<double(const double, const double)>&
        max_distance_from_linear,
    const std::function<bool(
        const double /*ti0*/,
        const double /*ti1*/,
        const double /*min_distance*/,
        const bool /*cm_no_zero_toi*/,
        double& /*toi*/)>& linear_ccd,
    double& toi,
    const double min_distance,
    const double tmax,
    const double conservative_rescaling)
{
    const double distance_t0 = distance(0);
    if (check_initial_distance(distance_t0, min_distance, toi)) {
        return true;
    }
    assert(distance_t0 > min_distance);

    double ti0 = 0;
    std::stack<double> ts;

// Initialize the stack of ts
#ifdef USE_FIXED_PIECES
    for (int i = FIXED_NUM_PIECES; i > 0; i--) {
        ts.push(i / double(FIXED_NUM_PIECES) * tmax);
    }
    int num_subdivisions = FIXED_NUM_PIECES;
#else
    ts.push(tmax);
    int num_subdivisions = 1;
#endif

    while (!ts.empty()) {
        const double ti1 = ts.top();

        const double distance_ti0 = distance(ti0);

        // If distance has decreased by a factor and the toi is not near zero,
        // then we can call this a collision.
        if (distance_ti0 < (1 - conservative_rescaling) * distance_t0
            && ti0 >= TightInclusionCCD::SMALL_TOI) {
            toi = ti0;
            logger().trace(
                "Distance small enough distance_ti0={:g}; toi={:g}",
                distance_ti0, toi);
            return true;
        }

        double min_distance_linear = max_distance_from_linear(ti0, ti1);

#ifndef USE_FIXED_PIECES
        // Check if the minimum distance is too large and we need to subdivide
        // (Large distances cause the slow CCD)
        if ((min_distance_linear
             >= std::min((1 - conservative_rescaling) * distance_ti0, 0.01))
            && (num_subdivisions < MAX_NUM_SUBDIVISIONS || ti0 == 0)) {
            logger().trace(
                "Subdividing at ti=[{:g}, {:g}] min_distance={:g} distance_ti0={:g}",
                ti0, ti1, min_distance_linear, distance_ti0);
            ts.push((ti1 + ti0) / 2);
            num_subdivisions++;
            continue;
        }
#endif

        min_distance_linear += min_distance;

        const bool is_impacting = linear_ccd(
            ti0, ti1, min_distance_linear, /*cm_no_zero_toi=*/ti0 == 0, toi);

        logger().trace(
            "Evaluated at ti=[{:g}, {:g}] min_distance={:g} distance_ti0={:g}; result={}{}",
            ti0, ti1, min_distance_linear, distance_ti0, is_impacting,
            is_impacting ? fmt::format(" toi={:g}", (ti1 - ti0) * toi + ti0)
                         : "");

        if (is_impacting) {
            toi = (ti1 - ti0) * toi + ti0;
            if (toi == 0) {
                // This is impossible because distance_t0 > min_distance
                ts.push((ti1 + ti0) / 2);
                num_subdivisions++;
                continue;
            }
            return true;
        }

        ts.pop();
        ti0 = ti1;
    }

    return false;
}

// ============================================================================

bool NonlinearCCD::point_point_ccd(
    const NonlinearTrajectory& p0,
    const NonlinearTrajectory& p1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    return conservative_piecewise_linear_ccd(
        [&](const double t) {
            return sqrt(point_point_distance(p0(t), p1(t)));
        },
        [&](const double t0, const double t1) {
            return p0.max_distance_from_linear(t0, t1)
                + p1.max_distance_from_linear(t0, t1);
        },
        [&](const double ti0, const double ti1, const double _min_distance,
            const bool cm_no_zero_toi, double& _toi) {
            return m_narrow_phase_ccd->point_point_ccd(
                p0(ti0), p1(ti0), p0(ti1), p1(ti1), _toi, _min_distance,
                /*tmax=*/1.0, cm_no_zero_toi);
        },
        toi, min_distance, tmax, conservative_rescaling);
}

bool NonlinearCCD::point_edge_ccd(
    const NonlinearTrajectory& p,
    const NonlinearTrajectory& e0,
    const NonlinearTrajectory& e1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    return conservative_piecewise_linear_ccd(
        [&](const double t) {
            return sqrt(point_edge_distance(p(t), e0(t), e1(t)));
        },
        [&](const double t0, const double t1) {
            return p.max_distance_from_linear(t0, t1)
                + std::max(
                       e0.max_distance_from_linear(t0, t1),
                       e1.max_distance_from_linear(t0, t1));
        },
        [&](const double ti0, const double ti1, const double _min_distance,
            const bool cm_no_zero_toi, double& _toi) {
            return m_narrow_phase_ccd->point_edge_ccd(
                p(ti0), e0(ti0), e1(ti0), p(ti1), e0(ti1), e1(ti1), _toi,
                _min_distance, /*tmax=*/1.0, cm_no_zero_toi);
        },
        toi, min_distance, tmax, conservative_rescaling);
}

bool NonlinearCCD::edge_edge_ccd(
    const NonlinearTrajectory& ea0,
    const NonlinearTrajectory& ea1,
    const NonlinearTrajectory& eb0,
    const NonlinearTrajectory& eb1,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    return conservative_piecewise_linear_ccd(
        [&](const double t) {
            return sqrt(edge_edge_distance(ea0(t), ea1(t), eb0(t), eb1(t)));
        },
        [&](const double t0, const double t1) {
            return std::max(
                       ea0.max_distance_from_linear(t0, t1),
                       ea1.max_distance_from_linear(t0, t1))
                + std::max(
                       eb0.max_distance_from_linear(t0, t1),
                       eb1.max_distance_from_linear(t0, t1));
        },
        [&](const double ti0, const double ti1, const double _min_distance,
            const bool cm_no_zero_toi, double& _toi) {
            return m_narrow_phase_ccd->edge_edge_ccd(
                ea0(ti0), ea1(ti0), eb0(ti0), eb1(ti0), //
                ea0(ti1), ea1(ti1), eb0(ti1), eb1(ti1), //
                _toi, _min_distance, /*tmax=*/1.0, cm_no_zero_toi);
        },
        toi, min_distance, tmax, conservative_rescaling);
}

bool NonlinearCCD::point_triangle_ccd(
    const NonlinearTrajectory& p,
    const NonlinearTrajectory& t0,
    const NonlinearTrajectory& t1,
    const NonlinearTrajectory& t2,
    double& toi,
    const double min_distance,
    const double tmax) const
{
    return conservative_piecewise_linear_ccd(
        [&](const double t) {
            return sqrt(point_triangle_distance(p(t), t0(t), t1(t), t2(t)));
        },
        [&](const double ti0, const double ti1) {
            return p.max_distance_from_linear(ti0, ti1)
                + std::max(
                       { t0.max_distance_from_linear(ti0, ti1),
                         t1.max_distance_from_linear(ti0, ti1),
                         t2.max_distance_from_linear(ti0, ti1) });
        },
        [&](const double ti0, const double ti1, const double _min_distance,
            const bool cm_no_zero_toi, double& _toi) {
            return m_narrow_phase_ccd->point_triangle_ccd(
                p(ti0), t0(ti0), t1(ti0), t2(ti0), //
                p(ti1), t0(ti1), t1(ti1), t2(ti1), //
                _toi, _min_distance, /*tmax=*/1.0, cm_no_zero_toi);
        },
        toi, min_distance, tmax, conservative_rescaling);
}

} // namespace ipc
