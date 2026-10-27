// PROTOTYPE ONLY: runtime switches for comparing edge-edge mollifiers.
#pragma once

#include <cmath>
#include <cstdlib>
#include <cstring>

namespace ipc {

/// Read once per process from the environment:
///   ESP_EE_MOLLIFIER=gcp|smoothstep|gcp-fixed|gcp-sqrt (default gcp)
///     gcp        x = c^2 / (eps d^2)       current (width sqrt(eps) d)
///     gcp-fixed  x = c^2 / l^2             same shape, fixed width l
///     gcp-sqrt   x = c^2 / (eps d L)       same shape, width sqrt(eps d L)
///     smoothstep H(c / l), H = 10t^3 - 15t^4 + 6t^5
///     smoothstep-sqrt  H(c / w), w^2 = eps d L
///     smoothstep-rel   H(c / w), w^2 = eps_rel d^2 (C2 shape, width
///                      relative to the gap; eps_rel = 1 by default, w = d)
///   with c^2 = d_PE^2 - d^2 the endpoint clearance.
///   ESP_EE_MOLLIFIER_EPS=<float>      (eps_rel for smoothstep-rel, default 1:
///                                      the EE sample fades in over the gap d,
///                                      the scale on which its endpoint vertex
///                                      loses potential; w <= d keeps the
///                                      barrier near endpoints)
///   ESP_EE_MOLLIFIER_ELL=<meters>     (fixed width l, default 150e-6)
///   ESP_EE_MOLLIFIER_L=<meters>       (gcp-sqrt length L, default 1.5e-3)
///   ESP_MOLLIFIER_ORDER=<int>         (override the per-barrier power; 0 =
///   keep) ESP_EE_AREA_REF=<meters>|dbar     (0 = off: scale each EE weight by
///   the
///                                      smooth min (1 + (a/rho)^2)^(-1/2), with
///                                      a = rest face area / 3, rho = r^2
///                                      (faces) a = rest rod length / 2, rho =
///                                      r (rods); "dbar" uses each
///                                      ESPParameters' dbar)
///   ESP_LOCAL_TRANSFER=1              (default for
///                                      ESPParameters::local_quadrature_transfer)
///   ESP_EE_WEIGHT=<float>             (EE sample weight factor, default 1.5:
///                                      1.5 * cubic_spline(0) = 1 vertex
///                                      sample)
struct EEMollifierConfig {
    bool smoothstep = false;
    bool fixed = false; // gcp-fixed
    bool sqrt_width = false; // gcp-sqrt, smoothstep-sqrt
    bool rel_width = false; // smoothstep-rel
    // smoothstep-rel width^2 = rel_eps * d^2 (ESP_EE_MOLLIFIER_EPS)
    double rel_eps = 1.0;
    double ell = 150e-6;
    double L = 1.5e-3;
    int order_override = 0;
    double ee_weight = 1.5;
    double area_ref = 0;
    bool area_ref_dbar = false; // ESP_EE_AREA_REF=dbar
    bool local_transfer_default = false;
};

/// Area-consistent EE weight factor for an element: smooth min(1, rho / a).
/// `a` is the rest area (faces: A/3) or length (rods: L/2) per vertex sample,
/// `rho` is r^2 (faces) or r (rods). Returns 1 when the feature is off.
inline double ee_area_scale_factor(double a, double r, bool face)
{
    if (r <= 0) {
        return 1.0;
    }
    const double t = a / (face ? r * r : r);
    return 1.0 / std::sqrt(1.0 + t * t);
}

inline const EEMollifierConfig& ee_mollifier_config()
{
    static const EEMollifierConfig config = [] {
        EEMollifierConfig c;
        if (const char* s = std::getenv("ESP_EE_MOLLIFIER")) {
            c.smoothstep = std::strcmp(s, "smoothstep") == 0
                || std::strcmp(s, "smoothstep-sqrt") == 0
                || std::strcmp(s, "smoothstep-rel") == 0;
            c.rel_width = std::strcmp(s, "smoothstep-rel") == 0;
            c.fixed = std::strcmp(s, "gcp-fixed") == 0;
            c.sqrt_width = std::strcmp(s, "gcp-sqrt") == 0
                || std::strcmp(s, "smoothstep-sqrt") == 0;
        }
        if (const char* s = std::getenv("ESP_EE_MOLLIFIER_EPS")) {
            c.rel_eps = std::atof(s);
        }
        if (const char* s = std::getenv("ESP_EE_MOLLIFIER_L")) {
            c.L = std::atof(s);
        }
        if (const char* s = std::getenv("ESP_EE_MOLLIFIER_ELL")) {
            c.ell = std::atof(s);
        }
        if (const char* s = std::getenv("ESP_EE_AREA_REF")) {
            c.area_ref_dbar = std::strcmp(s, "dbar") == 0;
            c.area_ref = c.area_ref_dbar ? 0.0 : std::atof(s);
        }
        if (const char* s = std::getenv("ESP_LOCAL_TRANSFER")) {
            c.local_transfer_default = std::strcmp(s, "1") == 0;
        }
        if (const char* s = std::getenv("ESP_EE_WEIGHT")) {
            c.ee_weight = std::atof(s);
        }
        if (const char* s = std::getenv("ESP_MOLLIFIER_ORDER")) {
            c.order_override = std::atoi(s);
        }
        return c;
    }();
    return config;
}

} // namespace ipc
