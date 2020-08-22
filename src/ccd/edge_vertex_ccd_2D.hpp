#pragma once

#include <Eigen/Core>

// namespace ipc {
namespace ccd {

/// Computes the time of impact between an \f$edge_{ij}\f$ and vertex
/// \f$k\f$.
///
/// @param[in]  V_{i,j,k} Vertex positions
/// @param[in]  U_{i,j,k} Vertex displacements
/// @param[out] toi       Time of FIRST impact
///
/// @return true if an impact happens at time \f$ t \in [0, 1]\f$
bool compute_edge_vertex_time_of_impact(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    double& toi,
    double& alpha);

/// @brief Given toi returns the position \f$\alpha\f$ along the
/// \f$edge_{ij}\f$ where the impact between vertex \f$k\f$ and
/// \f$edge_{ij}\f$  takes place.
///
/// @param[in] V_{i,j,k} Vertex positions
/// @param[in] U_{i,j,k} Vertex displacements
/// @param[in] toi       Time of impact
/// @param[out] alpha Position along the \f$edge_{ij}\f$ where the impact
/// takes place
///
///
bool temporal_parameterization_to_spatial(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    const double toi,
    double& alpha);

} // namespace ccd
// } // namespace ipc
