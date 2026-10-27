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
    /// Requires use_near_far=true and a face rule with positive weight at
    /// every corner (the corners donate the weight).
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

    /// Fixed face quadrature rule (3D), with weights summing to 1. The
    /// default is the vertex rule: the three corners with weight 1/3 each.
    /// Unused under set_auto_lattice(); see face_rule(). quad_order does not
    /// affect 3D.
    const FaceQuadRule& quad_rule() const { return m_face_quad_rule; }

    /// Weight of the rule's point at each corner of the face (0 if the rule
    /// has no point there). Corner points are integrated with the shared
    /// per-vertex samples, so each vertex is evaluated once, not once per
    /// incident face.
    const std::array<double, 3>& corner_weights() const
    {
        return m_corner_weights;
    }

    /// The rule's points that are not face corners (edge and interior
    /// points), each with its own per-face sample.
    const FaceQuadRule& non_corner_rule() const { return m_non_corner_rule; }

    /// @brief Set the face quadrature rule (3D).
    /// @param rule Quadrature points in barycentric coordinates with their
    ///        weights. The weights are normalized to sum to 1, so every rule
    ///        weighs a face's samples on the same scale as edge-edge samples.
    /// @throws std::invalid_argument if the rule is empty, if a weight is
    ///         negative or not finite, if the weights sum to zero, or if a
    ///         barycentric coordinate is negative (point outside the face) or
    ///         not finite. A negative weight makes the barrier energy unbounded
    ///         below as its point approaches contact.
    /// @note Logs a warning for each point whose barycentric coordinates do
    ///       not sum to 1 (within 1e-10).
    void set_quad_rule(FaceQuadRule rule)
    {
        if (rule.empty()) {
            log_and_throw_error<std::invalid_argument>(
                "Face quadrature rule is empty.");
        }
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
        if (weight_sum == 0) {
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
        m_corner_weights = { { 0, 0, 0 } };
        m_non_corner_rule.clear();
        for (auto& qp : rule) {
            qp.weight /= weight_sum;
            const int c = corner_index(qp.lambda);
            if (c >= 0) {
                m_corner_weights[c] += qp.weight;
            } else {
                m_non_corner_rule.push_back(qp);
            }
        }
        m_face_quad_rule = std::move(rule);
        m_lattice_spacing = 0; // a fixed rule replaces the auto lattice
        m_lattice_non_corner.clear();
    }

    /// @brief The vertex rule: the three face corners, weight 1/3 each.
    static FaceQuadRule vertex_rule()
    {
        return { { { { 1, 0, 0 } }, 1. / 3 },
                 { { { 0, 1, 0 } }, 1. / 3 },
                 { { { 0, 0, 1 } }, 1. / 3 } };
    }

    /// @brief Composite vertex rule on m x m sub-triangles.
    ///
    /// Splits the face into m^2 congruent sub-triangles and applies the
    /// vertex rule to each, so its points are the lattice (i, j, k) / m with
    /// i + j + k = m. A point's weight is the number of sub-triangles
    /// sharing it (1 at a corner, 3 on an edge, 6 inside) over 3 m^2. The
    /// weights are positive, the corners and edge points are included, and
    /// the point spacing is the edge length over m. Exact for linear
    /// functions only, but robust for the near-contact integrand, which is
    /// a barrier peak about a gap wide rather than a smooth polynomial.
    /// m = 1 is the vertex rule.
    static FaceQuadRule lattice_rule(const int m)
    {
        if (m < 1) {
            log_and_throw_error<std::invalid_argument>(
                "Lattice rule subdivisions {} must be at least 1.", m);
        }
        FaceQuadRule rule;
        const double inv = 1.0 / m, w = 1.0 / (3.0 * m * m);
        for (int i = 0; i <= m; i++) {
            for (int j = 0; i + j <= m; j++) {
                const int k = m - i - j;
                const int zeros = int(i == 0) + int(j == 0) + int(k == 0);
                const double count = zeros == 2 ? 1 : (zeros == 1 ? 3 : 6);
                rule.push_back(
                    { { { i * inv, j * inv, k * inv } }, count * w });
            }
        }
        return rule;
    }

    /// @brief Choose each face's lattice rule from its rest area (3D).
    ///
    /// Face f uses lattice_rule(m) with m = ceil(sqrt(4 A_f / sqrt(3)) /
    /// spacing), clamped to [1, max_subdivisions]: the side of an
    /// equilateral triangle of area A_f over the target point spacing. The
    /// rule depends on rest geometry only, so it is fixed per face and adds
    /// no non-smoothness in the deformed positions. It does change
    /// discontinuously with the rest shape, so use a fixed rule
    /// (set_quad_rule()) when differentiating with respect to the rest
    /// shape.
    /// @param spacing Target point spacing, e.g. a multiple of dhat.
    /// @param max_subdivisions Cap on m; a face uses at most
    ///        (m + 1)(m + 2) / 2 points.
    void set_auto_lattice(const double spacing, const int max_subdivisions = 16)
    {
        if (!std::isfinite(spacing) || spacing <= 0) {
            log_and_throw_error<std::invalid_argument>(
                "Lattice spacing {} must be positive and finite.", spacing);
        }
        if (max_subdivisions < 1) {
            log_and_throw_error<std::invalid_argument>(
                "Maximum lattice subdivisions {} must be at least 1.",
                max_subdivisions);
        }
        m_lattice_spacing = spacing;
        m_lattice_non_corner.assign(max_subdivisions + 1, {});
        for (int m = 2; m <= max_subdivisions; m++) {
            for (const auto& qp : lattice_rule(m)) {
                if (corner_index(qp.lambda) < 0) {
                    m_lattice_non_corner[m].push_back(qp);
                }
            }
        }
    }

    /// Target point spacing of the automatic lattice rules, or 0 if faces use
    /// the fixed rule (quad_rule()).
    double lattice_spacing() const { return m_lattice_spacing; }

    /// Cap on the automatic lattice subdivisions (0 without auto lattice).
    int max_lattice_subdivisions() const
    {
        return std::max(0, int(m_lattice_non_corner.size()) - 1);
    }

    /// Lattice subdivisions m for a face of rest area @p rest_area (auto
    /// lattice only).
    int lattice_subdivisions(const double rest_area) const
    {
        const double side = std::sqrt(4 * rest_area / std::sqrt(3.0));
        const int m = int(std::ceil(side / m_lattice_spacing - 1e-12));
        return std::clamp(m, 1, max_lattice_subdivisions());
    }

    /// A face's rule, split into corner weights and other points.
    struct FaceRule {
        std::array<double, 3> corner_weights;
        const FaceQuadRule* non_corner;
    };

    /// @brief The rule of a face of rest area @p rest_area: lattice_rule(m)
    /// with m from the area under set_auto_lattice(), else quad_rule().
    FaceRule face_rule(const double rest_area) const
    {
        if (m_lattice_spacing <= 0) {
            return { m_corner_weights, &m_non_corner_rule };
        }
        const int m = lattice_subdivisions(rest_area);
        const double c = 1.0 / (3.0 * m * m);
        return { { { c, c, c } }, &m_lattice_non_corner[m] };
    }

    /// Whether any face may have points other than its corners.
    bool has_non_corner_points() const
    {
        return !m_non_corner_rule.empty() || max_lattice_subdivisions() > 1;
    }

private:
    /// Index of the corner a point sits on, or -1.
    static int corner_index(const std::array<double, 3>& lambda)
    {
        constexpr double tol = 1e-12;
        for (int i = 0; i < 3; i++) {
            if (std::abs(lambda[i] - 1) <= tol
                && std::abs(lambda[(i + 1) % 3]) <= tol
                && std::abs(lambda[(i + 2) % 3]) <= tol) {
                return i;
            }
        }
        return -1;
    }

    const double m_dbar_factor;
    FaceQuadRule m_face_quad_rule = vertex_rule();
    std::array<double, 3> m_corner_weights { { 1. / 3, 1. / 3, 1. / 3 } };
    FaceQuadRule m_non_corner_rule;
    double m_lattice_spacing = 0;
    std::vector<FaceQuadRule> m_lattice_non_corner; ///< Indexed by m
};

} // namespace ipc
