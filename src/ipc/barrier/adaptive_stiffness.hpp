// Fucntions for computing the initial and updated barrier stiffnesses.

#pragma once

#include <ipc/collision_mesh.hpp>
#include <ipc/barrier/barrier.hpp>
#include <ipc/candidates/collision_stencil.hpp>

#include <Eigen/Core>

namespace ipc {

/// @brief Compute an inital barrier stiffness using the barrier potential gradient.
/// @param[in] bbox_diagonal Length of the diagonal of the bounding box of the scene.
/// @param[in] barrier Barrier function.
/// @param[in] dhat Activation distance of the barrier.
/// @param[in] average_mass Average mass of all bodies.
/// @param[in] grad_energy Gradient of the elasticity energy function.
/// @param[in] grad_barrier Gradient of the barrier potential.
/// @param[out] max_barrier_stiffness Maximum stiffness of the barrier.
/// @param[in] min_barrier_stiffness_scale Scale used to premultiply the minimum barrier stiffness.
/// @param[in] dmin Minimum distance between elements.
/// @return The initial barrier stiffness.
double initial_barrier_stiffness(
    const double bbox_diagonal,
    const Barrier& barrier,
    const double dhat,
    const double average_mass,
    Eigen::ConstRef<Eigen::VectorXd> grad_energy,
    Eigen::ConstRef<Eigen::VectorXd> grad_barrier,
    double& max_barrier_stiffness,
    const double min_barrier_stiffness_scale = 1e11,
    const double dmin = 0);

/// @brief Update the barrier stiffness if the distance is decreasing and less than dhat_epsilon_scale * diag.
/// @param[in] prev_min_distance Previous minimum distance between elements.
/// @param[in] min_distance Current minimum distance between elements.
/// @param[in] max_barrier_stiffness Maximum stiffness of the barrier.
/// @param[in] barrier_stiffness Current barrier stiffness.
/// @param[in] bbox_diagonal Length of the diagonal of the bounding box of the scene.
/// @param[in] dhat_epsilon_scale Update if distance is less than this fraction of the diagonal.
/// @param[in] dmin Minimum distance between elements.
/// @return The updated barrier stiffness.
double update_barrier_stiffness(
    const double prev_min_distance,
    const double min_distance,
    const double max_barrier_stiffness,
    const double barrier_stiffness,
    const double bbox_diagonal,
    const double dhat_epsilon_scale = 1e-9,
    const double dmin = 0);

/// @brief Compute the semi-implicit stiffness for a single collision.
/// @note See [Ando 2024] for details.
/// @param stencil Collision stencil.
/// @param vertices Vertex positions.
/// @param mass Vertex masses.
/// @param local_hess Local hessian of the elasticity energy function.
/// @param dmin Minimum distance between elements.
/// @return The semi-implicit stiffness.
/// @tparam StencilT The stencil type; any candidate or collision.
template <typename StencilT>
double semi_implicit_stiffness(
    const StencilT& stencil,
    Eigen::ConstRef<VectorMax12d> vertices,
    Eigen::ConstRef<VectorMax4d> mass,
    Eigen::ConstRef<MatrixMax12d> local_hess,
    const double dmin)
{
    const unsigned N = stencil.num_vertices();
    assert(vertices.size() % N == 0);
    const unsigned dim = stencil.dim(vertices.size());

    const VectorMax4d value = stencil.compute_coefficients(vertices);

    // Compute the contact normal (i.e., the vector from the )
    VectorMax3d normal = VectorMax3d::Zero(dim);
    for (unsigned i = 0; i < N; ++i) {
        normal += value[i] * vertices.segment(dim * i, dim);
    }

    // d²
    const double distance = normal.norm() - dmin;
    const double distance_sqr = distance * distance;

    // average mass: mᵢ = cᵀMc / ‖c‖²
    const double avg_mass =
        value.dot(mass.asDiagonal() * value) / value.squaredNorm();

    VectorMax12d w = VectorMax12d::Zero(dim * N);
    for (unsigned i = 0; i < N; ++i) {
        w.segment(dim * i, dim) = value[i] * normal;
    }
    w.normalize();

    return avg_mass / distance_sqr + w.dot(local_hess * w);
}

/// @brief Compute the semi-implicit stiffness's for all collisions.
/// @note See [Ando 2024] for details.
/// @param mesh Collision mesh.
/// @param vertices Vertex positions.
/// @param collisions Normal collisions or collision candidates.
/// @param vertex_masses Lumped vertex masses.
/// @param hess Hessian of the elasticity energy function.
/// @param dmin Minimum distance between elements.
/// @return The semi-implicit stiffness's.
template <typename StencilsT>
Eigen::VectorXd semi_implicit_stiffness(
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> vertices,
    const StencilsT& collisions,
    Eigen::ConstRef<Eigen::VectorXd> vertex_masses,
    const Eigen::SparseMatrix<double>& hess,
    const double dmin);

} // namespace ipc
