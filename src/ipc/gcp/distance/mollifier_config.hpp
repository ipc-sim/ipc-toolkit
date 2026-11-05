#pragma once

#include <cmath>
#include <cstdint>

namespace ipc {

/// @brief Shape of the 3D edge-edge mollifier.
///
/// An edge-edge sample fades in as the closest points move off the edge
/// endpoints. Each endpoint contributes a factor of its clearance
/// c² = d_PE² − d², where d_PE is the endpoint's distance to the other edge
/// and d the edge-edge distance.
enum class EEMollifierType : std::uint8_t {
    /// Smoothstep H(c / w), H(t) = 10t³ − 15t⁴ + 6t⁵ (C²), with a width
    /// relative to the gap, w² = rel_eps d². With rel_eps = 1 (w = d) the
    /// sample fades in over the same scale on which its endpoint vertex
    /// loses potential, the widest width that leaves no hole near the
    /// endpoints.
    SMOOTHSTEP,
    /// Original GCP mollifier: x = c² / (MOLLIFIER_THRESHOLD_EPS d²).
    GCP,
};

/// @brief Options of the 3D edge-edge mollifier (see EEMollifierType).
struct EEMollifierConfig {
    EEMollifierType type = EEMollifierType::SMOOTHSTEP;
    /// SMOOTHSTEP: squared width over squared gap, w² = rel_eps d².
    double rel_eps = 1.0;

    /// The original GCP mollifier, used by the GCP formulation.
    static EEMollifierConfig gcp()
    {
        EEMollifierConfig c;
        c.type = EEMollifierType::GCP;
        return c;
    }

    bool smoothstep() const { return type == EEMollifierType::SMOOTHSTEP; }
};

/// @brief Area-consistent edge-edge weight factor of an element:
/// (1 + (a / ρ)²)^(-1/2), a smooth min(1, ρ / a).
/// @param a Rest measure compared with the patch (see
///        ESPParameters::ee_area_scaling).
/// @param r Reference length; ρ = r² for faces and r for rods.
/// @param face Whether the element is a face.
/// @return The factor, or 1 when r <= 0.
inline double ee_area_scale_factor(double a, double r, bool face)
{
    if (r <= 0) {
        return 1.0;
    }
    const double t = a / (face ? r * r : r);
    return 1.0 / std::sqrt(1.0 + t * t);
}

} // namespace ipc
