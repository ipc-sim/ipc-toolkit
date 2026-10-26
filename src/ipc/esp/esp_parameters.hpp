#pragma once

#include <ipc/barrier/barrier.hpp>
#include <ipc/gcp/distance/mollifier_config.hpp>
#include <ipc/utils/logger.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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

    /// @param _dhat Activation distance of the barrier; must be positive and
    ///        finite.
    /// @param _dbar_factor Near/far split parameter dbar / dhat, in
    ///        [0, 2]. The far part of the barrier vanishes identically at 2,
    ///        so the split tends continuously to the unsplit potential; the
    ///        edge-edge weights use the support min(dbar, dhat) (ee_support).
    ///        0 disables edge-edge contact.
    /// @param _quad_order Gauss-Lobatto edge quadrature order (2D); see
    ///        quad_order.
    /// @throws std::invalid_argument if _dhat is not positive and finite, if
    ///         _dbar_factor is not in [0, 2], or if _quad_order is negative.
    ESPParameters(
        const double _dhat,
        const double _dbar_factor = 0.2,
        const int _quad_order = 1,
        bool _area_weights = true,
        const IntegrationType _integration_type = IntegrationType::NORMAL)
        : dhat(_dhat)
        , dbar(_dbar_factor * dhat)
        , quad_order(_quad_order)
        , area_weights(_area_weights)
        , integration_type(_integration_type)
        , m_dbar_factor(_dbar_factor)
    {
        // A zero dhat silently disables contact, and a negative or non-finite
        // one produces negative or NaN broad-phase inflation radii.
        if (!std::isfinite(_dhat) || _dhat <= 0) {
            log_and_throw_error<std::invalid_argument>(
                "dhat {} must be positive and finite.", _dhat);
        }
        if (_quad_order < 0) {
            log_and_throw_error<std::invalid_argument>(
                "quad_order {} must be non-negative.", _quad_order);
        }
        if (std::isnan(_dbar_factor) || _dbar_factor < 0 || _dbar_factor > 2) {
            log_and_throw_error<std::invalid_argument>(
                "dbar_factor {} is not in [0, 2].", _dbar_factor);
        } else if (_dbar_factor == 0) {
            logger().warn("dbar_factor = 0 disables edge-edge contact.");
        }
    }

    const double dhat;
    const double dbar;

    /// Barrier evaluated on each pair's unsquared distance d with activation
    /// distance dhat, in 2D and 3D. With an InversePowerBarrier of power p, the
    /// 3D edge-edge mollifier weights are raised to the power max(1, round(p) +
    /// 1).
    std::shared_ptr<Barrier> barrier =
        std::make_shared<NormalizedClampedLogBarrier<>>();
    /// Gauss-Lobatto edge quadrature order (2D): order n uses n + 1 points and
    /// must be at least 1 in 2D (ESPCollisions::build() checks this). Unused
    /// in 3D, where any non-negative value is accepted.
    const int quad_order;
    /// Weight each face by its area (3D), or each edge by its area
    /// (CollisionMesh::edge_area(), 2D). In 3D, an edge in no face is weighted
    /// by its edge area and an isolated vertex by its vertex area. If false,
    /// every weight is 1.
    bool area_weights;
    /// Experimental 3D near-field quadrature: transfer weight locally from
    /// edge endpoints to EE samples instead of renormalizing all samples.
    /// Requires use_near_far=true and the default vertex quadrature rule.
    /// Currently supports frictionless contact only.
    /// Has no effect in 2D or when dbar_factor=0. Disabled by default.
    bool local_quadrature_transfer =
        ee_mollifier_config().local_transfer_default; // PROTOTYPE env default
    const IntegrationType integration_type;

    /// @brief Near/far split parameter dbar / dhat, in [0, 2].
    double dbar_factor() const { return m_dbar_factor; }

    /// @brief Support of the edge-edge (edge point) weights: min(dbar, dhat).
    /// Edge-edge pairs farther apart than this get no quadrature point.
    double ee_support() const { return std::min(dbar, dhat); }

    /// Face quadrature rule (3D). Empty (default) integrates over the face
    /// vertices instead; quad_order does not affect 3D.
    const FaceQuadRule& quad_rule() const { return m_face_quad_rule; }

    /// @brief Set the face quadrature rule (3D).
    /// @param rule Quadrature points in barycentric coordinates with their
    ///        weights. Empty integrates over the face vertices instead.
    /// @throws std::invalid_argument if a weight is negative or not finite, if
    ///         the weights sum to zero, or if a barycentric coordinate is
    ///         negative (point outside the face) or not finite. A negative
    ///         weight makes the barrier energy unbounded below as its point
    ///         approaches contact.
    /// @note Logs a warning for each point whose barycentric coordinates do
    ///       not sum to 1 (within 1e-10).
    void set_quad_rule(FaceQuadRule rule)
    {
        double weight_sum = 0;
        std::vector<size_t> not_affine; // points whose coordinates sum != 1
        for (size_t i = 0; i < rule.size(); i++) {
            const double weight = rule[i].weight;
            if (!std::isfinite(weight) || weight < 0) {
                log_and_throw_error<std::invalid_argument>(
                    "Face quadrature point {} has negative or non-finite "
                    "weight {}.",
                    i, weight);
            }
            weight_sum += weight;

            double lambda_sum = 0;
            for (const double lambda : rule[i].lambda) {
                if (!std::isfinite(lambda) || lambda < 0) {
                    log_and_throw_error<std::invalid_argument>(
                        "Face quadrature point {} has negative or "
                        "non-finite barycentric coordinate {}.",
                        i, lambda);
                }
                lambda_sum += lambda;
            }
            if (std::abs(lambda_sum - 1) > 1e-10) {
                not_affine.push_back(i);
            }
        }
        if (!rule.empty() && weight_sum == 0) {
            log_and_throw_error<std::invalid_argument>(
                "Face quadrature weights sum to zero.");
        }
        for (const size_t i : not_affine) {
            const auto& l = rule[i].lambda;
            logger().warn(
                "Barycentric coordinates of face quadrature point {} sum to "
                "{}, not 1.",
                i, l[0] + l[1] + l[2]);
        }
        m_face_quad_rule = std::move(rule);
    }

private:
    const double m_dbar_factor;
    FaceQuadRule m_face_quad_rule;
};

} // namespace ipc
