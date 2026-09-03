#pragma once

#include "esp_parameters.hpp"

#include <ipc/collision_mesh.hpp>

#include <Eigen/Core>

namespace ipc {

/// Manages per-vertex dhat values for adaptive barrier support sizing.
/// dhat is defined only at vertices; edge/face values are linearly
/// interpolated.
class AdaptiveSupport {
public:
    /// Construct from rest mesh and positions. Computes and stores per-vertex
    /// dhat values.
    AdaptiveSupport(
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
        const ESPParameters& params);

    /// Get dhat value at a vertex.
    double vertex(index_t vertex_id) const;

    /// Get interpolated dhat on edge at barycentric parameter t in [0,1].
    /// t=0 corresponds to edge vertex 0, t=1 to edge vertex 1.
    double edge(index_t edge_id, double t) const;

    /// Get interpolated dhat on face at barycentric coordinates (u, v).
    /// u, v are barycentric coords; third coord is 1-u-v.
    double face(index_t face_id, double u, double v) const;

    /// Scale all per-vertex dhat values by a factor.
    void scale(double factor) { m_values *= factor; }

    /// Multiplicative reduction factor applied to primitive vertex dhat values
    /// each time they are found in a non-zero collision pair.
    double zeta = 0.8;

private:
    Eigen::VectorXd m_values;
    const CollisionMesh* m_mesh;
};

} // namespace ipc
