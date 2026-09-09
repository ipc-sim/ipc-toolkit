#pragma once

#include <ipc/barrier/barrier.hpp>
#include <ipc/utils/logger.hpp>

#include <array>
#include <atomic>
#include <cstdint>
#include <limits>
#include <memory>

namespace ipc {

/// A single face quadrature point in barycentric coordinates with its weight.
struct FaceQuadPoint {
    std::array<double, 3> lambda; ///< Barycentric coordinates (sum = 1)
    double weight;
};
using FaceQuadRule = std::vector<FaceQuadPoint>;

struct ESPParameters {
    enum class IntegrationType : std::uint8_t {
        BRUTE_FORCE, ///< Integrate all pairs with no obstacle filtering
        NORMAL, ///< Filter obstacle-obstacle pairs; skip primitives with only
                ///< obstacle candidates
        NO_OBST ///< Skip obstacle sources entirely, may miss collisions!
    };

    ESPParameters(
        const double _dhat,
        const double dbar_factor_value = 1.0,
        const int _quad_order = 1,
        bool _area_weights = true,
        const IntegrationType _integration_type = IntegrationType::NORMAL)
        : dhat(_dhat)
        , dbar(dbar_factor_value * dhat)
        , dbar_factor_value(dbar_factor_value)
        , quad_order(_quad_order)
        , area_weights(_area_weights)
        , integration_type(_integration_type)
    {
        if (quad_order > 14) {
            throw std::invalid_argument(
                "Quadrature order " + std::to_string(quad_order)
                + ">14 is too large.");
        } else if (quad_order == 6 || quad_order == 8) {
            logger().error(
                "Quadrature orders 6 and 8 has negative vertex weights.");
        } else if (quad_order >= 10 && quad_order <= 12) {
            logger().warn(
                "Quadrature orders 10-12 are not implemented, and instead use order 13.");
        } else if (quad_order == 1) {
            logger().warn(
                "Quadrature order 1 is equivalent to vertex quadrature.");
        }
    }

    const double dhat;
    const double dbar;
    const double dbar_factor_value;

    /// Barrier function used in 3D collision evaluation.
    std::shared_ptr<Barrier> barrier =
        std::make_shared<NormalizedClampedLogBarrier<>>();
    const int quad_order;
    bool area_weights;
    const IntegrationType integration_type;

    double dbar_factor() const { return dbar_factor_value; }

    const FaceQuadRule& get_quad_rule() const { return face_quad_rule; }

    /// Face quadrature rule. Empty (default) skips face quadrature entirely,
    /// matching the behaviour of quad_order == 0 in the old interface.
    FaceQuadRule face_quad_rule;

    /// Record a distance passed to the barrier; tracks the running minimum
    /// across all threads. Copies of this struct share the same tracker
    /// (shared_ptr), so pass-by-value sites still update the original.
    void record_dist(double d) const
    {
        auto& a = *m_min_dist_seen;
        double cur = a.load(std::memory_order_relaxed);
        while (d < cur
               && !a.compare_exchange_weak(cur, d, std::memory_order_relaxed)) {
        }
    }

    double min_dist_seen() const
    {
        return m_min_dist_seen->load(std::memory_order_relaxed);
    }

    void reset_min_dist() const
    {
        m_min_dist_seen->store(
            std::numeric_limits<double>::infinity(), std::memory_order_relaxed);
    }

private:
    std::shared_ptr<std::atomic<double>> m_min_dist_seen =
        std::make_shared<std::atomic<double>>(
            std::numeric_limits<double>::infinity());
};

} // namespace ipc
