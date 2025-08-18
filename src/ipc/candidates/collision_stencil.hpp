#pragma once

#include <ipc/config.hpp>
#include <ipc/ccd/default_narrow_phase_ccd.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <array>
#include <limits>
#include <ostream>

namespace ipc {

/// @brief A stencil representing a collision between at most four vertices.
class CollisionStencil {
public:
    virtual ~CollisionStencil() = default;

    /// @brief Get the number of vertices in the collision stencil.
    virtual int num_vertices() const = 0;

    /// @brief Get the dimension of the collision stencil.
    /// @param ndof Number of degrees of freedom in the stencil.
    /// @return The dimension of the collision stencil.
    int dim(const int ndof) const
    {
        assert(ndof % num_vertices() == 0);
        return ndof / num_vertices();
    }

    /// @brief Get the vertex IDs of the collision stencil.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return The vertex IDs of the collision stencil. Size is always 4, but elements i > num_vertices() are -1.
    virtual std::array<index_t, 4> vertex_ids(
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const = 0;

    /// @brief Get the vertex attributes of the collision stencil.
    /// @param vertices Vertex attributes
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return The vertex positions of the collision stencil. Size is always 4, but elements i > num_vertices() are NaN.
    std::array<VectorMax3d, 4> vertices(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        constexpr double NaN = std::numeric_limits<double>::signaling_NaN();

        const auto vertex_ids = this->vertex_ids(edges, faces);

        std::array<VectorMax3d, 4> stencil_vertices;
        for (int i = 0; i < 4; i++) {
            if (vertex_ids[i] >= 0) {
                stencil_vertices[i] = vertices.row(vertex_ids[i]);
            } else {
                stencil_vertices[i].setConstant(vertices.cols(), NaN);
            }
        }

        return stencil_vertices;
    }

    /// @brief Select this stencil's DOF from the full matrix of DOF.
    /// @param X Full matrix of DOF (rowwise).
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return This stencil's DOF.
    VectorMax12d
    dof(Eigen::ConstRef<Eigen::MatrixXd> X,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        const int dim = X.cols();
        VectorMax12d x(num_vertices() * dim);
        const auto idx = vertex_ids(edges, faces);
        for (int i = 0; i < num_vertices(); i++) {
            x.segment(i * dim, dim) = X.row(idx[i]);
        }
        return x;
    }

    /// @brief Compute the distance of the stencil.
    /// @param vertices Collision mesh vertices
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return Distance of the stencil.
    double compute_distance(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return compute_distance(dof(vertices, edges, faces));
    }

    /// @brief Compute the distance gradient of the stencil w.r.t. the stencil's vertex positions.
    /// @param vertices Collision mesh vertices
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return Distance gradient of the stencil w.r.t. the stencil's vertex positions.
    VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return compute_distance_gradient(dof(vertices, edges, faces));
    }

    /// @brief Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    /// @param vertices Collision mesh vertices
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    MatrixMax12d compute_distance_hessian(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return compute_distance_hessian(dof(vertices, edges, faces));
    }

    /// @brief Compute the coefficients of the stencil s.t. \f$d(x) = \|\sum c_i \mathbf{x}_i\|^2\f$.
    /// @param vertices Collision mesh vertices
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return Coefficients of the stencil.
    VectorMax4d compute_coefficients(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return compute_coefficients(dof(vertices, edges, faces));
    }

    /// @brief Compute the normal of the stencil.
    /// @param vertices Collision mesh vertices
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param flip_if_negative If true, flip the normal if the point is on the negative side.
    /// @param sign If not nullptr, set to the sign of the normal before any flipping.
    /// @return Normal of the stencil.
    VectorMax3d compute_normal(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const bool flip_if_negative = true,
        double* sign = nullptr) const
    {
        return compute_normal(
            dof(vertices, edges, faces), flip_if_negative, sign);
    }

    /// @brief Compute the Jacobian of the normal of the stencil.
    /// @param vertices Collision mesh vertices
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param flip_if_negative If true, flip the normal if the point is on the negative side.
    /// @return Jacobian of the normal of the stencil.
    MatrixMax<double, 3, 12> compute_normal_jacobian(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces,
        const bool flip_if_negative = true) const
    {
        return compute_normal_jacobian(
            dof(vertices, edges, faces), flip_if_negative);
    }

    // ----------------------------------------------------------------------
    // NOTE: The following functions take stencil vertices as output by dof()
    // ----------------------------------------------------------------------

    /// @brief Compute the distance of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @note positions can be computed as stencil.dof(vertices, edges, faces)
    /// @return Distance of the stencil.
    virtual double
    compute_distance(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the distance gradient of the stencil w.r.t. the stencil's vertex positions.
    /// @param positions Stencil's vertex positions.
    /// @note positions can be computed as stencil.dof(vertices, edges, faces)
    /// @return Distance gradient of the stencil w.r.t. the stencil's vertex positions.
    virtual VectorMax12d compute_distance_gradient(
        Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    /// @param positions Stencil's vertex positions.
    /// @note positions can be computed as stencil.dof(vertices, edges, faces)
    /// @return Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    virtual MatrixMax12d
    compute_distance_hessian(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the coefficients of the stencil s.t. d(x) = ‖∑ cᵢ xᵢ‖².
    /// @param positions Stencil's vertex positions.
    /// @return Coefficients of the stencil.
    virtual VectorMax4d
    compute_coefficients(Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the normal of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @param flip_if_negative If true, flip the normal if the point is on the negative side.
    /// @param sign If not nullptr, set to the sign of the normal before any flipping.
    /// @return Normal of the stencil.
    VectorMax3d compute_normal(
        Eigen::ConstRef<VectorMax12d> positions,
        bool flip_if_negative = true,
        double* sign = nullptr) const;

    /// @brief Compute the Jacobian of the normal of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @param flip_if_negative If true, flip the normal if the point is on the negative side.
    /// @return Jacobian of the normal of the stencil.
    MatrixMax<double, 3, 12> compute_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions,
        bool flip_if_negative = true) const;

    /// @brief Perform narrow-phase CCD on the candidate.
    /// @param[in] vertices_t0 Stencil vertices at the start of the time step.
    /// @param[in] vertices_t1 Stencil vertices at the end of the time step.
    /// @param[out] toi Computed time of impact (normalized).
    /// @param[in] min_distance Minimum separation distance between primitives.
    /// @param[in] tmax Maximum time (normalized) to look for collisions.
    /// @param[in] narrow_phase_ccd The narrow phase CCD algorithm to use.
    /// @return If the candidate had a collision over the time interval.
    virtual bool
    ccd(Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1,
        double& toi,
        const double min_distance = 0.0,
        const double tmax = 1.0,
        const NarrowPhaseCCD& narrow_phase_ccd =
            DEFAULT_NARROW_PHASE_CCD) const = 0;

    /// @brief Write the CCD query to a stream.
    /// @param out Stream to write to.
    /// @param vertices_t0 Stencil vertices at the start of the time step.
    /// @param vertices_t1 Stencil vertices at the end of the time step.
    /// @return The stream.
    std::ostream& write_ccd_query(
        std::ostream& out,
        Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1) const;

protected:
    /// @brief Compute the unnormalized normal of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @return Unnormalized normal of the stencil.
    virtual VectorMax3d compute_unnormalized_normal(
        Eigen::ConstRef<VectorMax12d> positions) const = 0;

    /// @brief Compute the Jacobian of the unnormalized normal of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @return Jacobian of the unnormalized normal of the stencil.
    virtual MatrixMax<double, 3, 12> compute_unnormalized_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const = 0;
};

} // namespace ipc