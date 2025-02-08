// Functions for computing the initial and updated barrier stiffnesses.

#include "adaptive_stiffness.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/candidates/candidates.hpp>
#include <ipc/collisions/normal/normal_collisions.hpp>

#include <algorithm> // std::min/max
#include <cassert>

namespace ipc {

double initial_barrier_stiffness(
    const double bbox_diagonal,
    const Barrier& barrier,
    const double dhat,
    const double average_mass,
    const Eigen::VectorXd& grad_energy,
    const Eigen::VectorXd& grad_barrier,
    double& max_barrier_stiffness,
    const double min_barrier_stiffness_scale,
    const double dmin)
{
    assert(average_mass > 0 && min_barrier_stiffness_scale > 0);
    assert(bbox_diagonal > 0);

    double dhat_squared = dhat * dhat;
    double dmin_squared = dmin * dmin;

    // Find a good initial value for κ
    double d0 = 1e-8 * bbox_diagonal + dmin;
    d0 *= d0;
    if (d0 - dmin_squared >= 2 * dmin * dhat + dhat_squared) {
        d0 = dmin * dhat + 0.5 * dhat_squared; // NOTE: this is untested
    }
    double min_barrier_stiffness = 4 * d0
        * barrier.second_derivative(
            d0 - dmin_squared, 2 * dmin * dhat + dhat_squared);
    min_barrier_stiffness =
        min_barrier_stiffness_scale * average_mass / min_barrier_stiffness;
    assert(std::isfinite(min_barrier_stiffness));

    max_barrier_stiffness = 100 * min_barrier_stiffness;

    double kappa = 1.0;
    if (grad_barrier.squaredNorm() > 0) {
        // If this value is negative it will be clamped to κ_min anyways
        kappa = -grad_barrier.dot(grad_energy) / grad_barrier.squaredNorm();
        assert(std::isfinite(kappa));
    }

    return std::clamp(kappa, min_barrier_stiffness, max_barrier_stiffness);
}

// Adaptive κ
double update_barrier_stiffness(
    const double prev_min_distance,
    const double min_distance,
    const double max_barrier_stiffness,
    const double barrier_stiffness,
    const double bbox_diagonal,
    const double dhat_epsilon_scale,
    const double dmin)
{
    // Is the barrier having a difficulty pushing the bodies apart?
    double dhat_epsilon = dhat_epsilon_scale * (bbox_diagonal + dmin);
    dhat_epsilon *= dhat_epsilon;
    if (prev_min_distance < dhat_epsilon && min_distance < dhat_epsilon
        && min_distance < prev_min_distance) {
        // Then increase the barrier stiffness.
        return std::min(max_barrier_stiffness, 2 * barrier_stiffness);
    }
    return barrier_stiffness;
}

// -----------------------------------------------------------------------------
//
// Based on `compute_stiffness()` in `src/cpp/barrier/barrier.cu` of
// (ppf-contact-solver)[https://github.com/st-tech/ppf-contact-solver]
//
// Original license:
// File: barrier.cu
// Author: Ryoichi Ando (ryoichi.ando@zozo.com)
// License: Apache v2.0
//

double semi_implicit_stiffness(
    const CollisionStencil& stencil,
    const std::array<long, 4>& vertex_ids,
    const VectorMax12d& vertices,
    const VectorMax4d& mass,
    const MatrixMax12d& local_hess,
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

template <typename StencilsT>
Eigen::VectorXd semi_implicit_stiffness(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const StencilsT& collisions,
    const Eigen::VectorXd& vertex_masses,
    const Eigen::SparseMatrix<double>& hess,
    const double dmin)
{
    const int dim = mesh.dim();
    assert(vertices.cols() == dim);     // Vertex positions must be 3D
    assert(hess.rows() == hess.cols()); // Hessian must be square
    // Hess and vertex_masses must have the same number of rows
    assert(hess.rows() == vertex_masses.size() * dim);
    // Hess can be either for the reduced or full mesh
    assert(hess.rows() == mesh.ndof() || hess.rows() == mesh.full_ndof());

    Eigen::VectorXd stiffnesses(collisions.size());

    for (size_t ci = 0; ci < collisions.size(); ci++) {
        const CollisionStencil& collision = collisions[ci];
        const unsigned N = collision.num_vertices();

        const VectorMax12d positions =
            collision.dof(vertices, mesh.edges(), mesh.faces());

        std::array<long, 4> vertex_ids =
            collision.vertex_ids(mesh.edges(), mesh.faces());
        if (hess.rows() == mesh.full_ndof()) {
            for (int i = 0; i < N; i++) {
                vertex_ids[i] = mesh.to_full_vertex_id(vertex_ids[i]);
            }
        }

        VectorMax4d local_mass(collision.num_vertices());
        for (unsigned i = 0; i < collision.num_vertices(); i++) {
            local_mass[i] = vertex_masses[vertex_ids[i]];
        }

        MatrixMax12d local_hess = MatrixMax12d::Zero(dim * N, dim * N);
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < N; ++j) {
                for (unsigned k = 0; k < dim; ++k) {
                    for (unsigned l = 0; l < dim; ++l) {
                        // NOTE: Assumes DOF are flattened in row-major order
                        local_hess(dim * i + k, dim * j + l) = hess.coeff(
                            dim * vertex_ids[i] + k, dim * vertex_ids[j] + l);
                    }
                }
            }
        }

        stiffnesses[ci] = semi_implicit_stiffness(
            collision, vertex_ids, positions, local_mass, local_hess, dmin);
    }

    return stiffnesses;
}

template Eigen::VectorXd semi_implicit_stiffness<NormalCollisions>(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const NormalCollisions& collisions,
    const Eigen::VectorXd& vertex_masses,
    const Eigen::SparseMatrix<double>& hess,
    const double dmin);

template Eigen::VectorXd semi_implicit_stiffness<Candidates>(
    const CollisionMesh& mesh,
    const Eigen::MatrixXd& vertices,
    const Candidates& collisions,
    const Eigen::VectorXd& vertex_masses,
    const Eigen::SparseMatrix<double>& hess,
    const double dmin);

// -----------------------------------------------------------------------------

} // namespace ipc
