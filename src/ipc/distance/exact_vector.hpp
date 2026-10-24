#pragma once

// This header depends on geogram's MultiPrecision PSM, a private dependency of
// ipc_toolkit, so only our own sources and tests include it. Our copy of the
// PSM lives in namespace ipc_geogram instead of GEO, so it never collides with
// a full geogram linked by a downstream project.
#include <MultiPrecision_psm.h>

namespace ipc {

/// @brief Exact scalar type (an adaptive-precision floating-point expansion).
using ExactReal = ipc_geogram::expansion_nt;

/// @brief Minimal exact 3D vector for sign predicates.
///
/// The MultiPrecision PSM only ships the scalar type, not geogram's vec3E, so
/// we provide the handful of operations the exact distance-type predicates
/// need. Every operation is exact, so the sign of any result is exact too.
struct ExactVector3 {
    ExactReal x, y, z;

    /// @brief Exact copy of a 2D or 3D point (z = 0 in 2D).
    static ExactVector3 from(const double* v, const int dim)
    {
        return { ExactReal(v[0]), ExactReal(v[1]),
                 ExactReal(dim == 3 ? v[2] : 0.0) };
    }

    ExactVector3 operator-(const ExactVector3& o) const
    {
        return { x - o.x, y - o.y, z - o.z };
    }

    ExactVector3 operator-() const { return { -x, -y, -z }; }

    ExactReal dot(const ExactVector3& o) const
    {
        return x * o.x + y * o.y + z * o.z;
    }

    ExactVector3 cross(const ExactVector3& o) const
    {
        return { y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x };
    }

    ExactReal squared_norm() const { return dot(*this); }
};

/// @brief Initialize the expansion arithmetic once (thread-safe).
///
/// geogram's MultiPrecision PSM computes its splitter and machine epsilon at
/// run time, so we must call this before the first exact computation.
inline void init_exact_arithmetic()
{
    static const bool initialized = [] {
        ipc_geogram::expansion::initialize();
        return true;
    }();
    (void)initialized;
}

} // namespace ipc
