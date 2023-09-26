#pragma once

#include <ipc/utils/eigen_ext.hpp>

#include <array>
#include <limits>

namespace ipc {

/// @brief A stencil representing a collision between at most four vertices.
class CollisionStencil {
public:
    virtual ~CollisionStencil() = default;

    /// @brief Get the number of vertices in the contact stencil.
    virtual int num_vertices() const = 0;

    /// @brief Get the vertex IDs of the contact stencil.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return The vertex IDs of the contact stencil. Size is always 4, but elements i > num_vertices() are -1.
    virtual std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const = 0;

    /// @brief Get the vertex attributes of the contact stencil.
    /// @tparam T Type of the attributes
    /// @param vertices Vertex attributes
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return The vertex positions of the contact stencil. Size is always 4, but elements i > num_vertices() are NaN.
    template <typename T>
    std::array<VectorMax3<T>, 4> vertices(
        const MatrixX<T>& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const
    {
        constexpr double NaN = std::numeric_limits<double>::signaling_NaN();

        const std::array<long, 4> vertex_ids = this->vertex_ids(edges, faces);

        std::array<VectorMax3<T>, 4> stencil_vertices;
        for (int i = 0; i < 4; i++) {
            if (vertex_ids[i] >= 0) {
                stencil_vertices[i] = vertices.row(vertex_ids[i]);
            } else {
                stencil_vertices[i].setConstant(vertices.cols(), T(NaN));
            }
        }

        return stencil_vertices;
    }

    /// @brief Select this stencil's DOF from the full matrix of DOF.
    /// @tparam T Type of the DOF
    /// @param X Full matrix of DOF (rowwise).
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return This constraint's DOF.
    template <typename T>
    VectorMax12<T>
    dof(const MatrixX<T>& X,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const
    {
        const int dim = X.cols();
        VectorMax12<T> x(num_vertices() * dim);
        const std::array<long, 4> idx = vertex_ids(edges, faces);
        for (int i = 0; i < num_vertices(); i++) {
            x.segment(i * dim, dim) = X.row(idx[i]);
        }
        return x;
    }

    /// @brief Compute the distance of the stencil.
    /// @param vertices Collision mesh vertex positions.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return Distance of the stencil.
    double compute_distance(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const;

    /// @brief Compute the distance gradient of the stencil w.r.t. the stencil's vertex positions.
    /// @param vertices Collision mesh vertex positions.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return Distance gradient of the stencil w.r.t. the stencil's vertex positions.
    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const;

    /// @brief Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    /// @param vertices Collision mesh vertex positions.
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const;

protected:
    /// @brief Compute the distance of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @return Distance of the stencil.
    virtual double compute_distance(const VectorMax12d& positions) const = 0;

    /// @brief Compute the distance gradient of the stencil w.r.t. the stencil's vertex positions.
    /// @param positions Stencil's vertex positions.
    /// @return Distance gradient of the stencil w.r.t. the stencil's vertex positions.
    virtual VectorMax12d
    compute_distance_gradient(const VectorMax12d& positions) const = 0;

    /// @brief Compute the distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    /// @param positions Stencil's vertex positions.
    /// @return Distance Hessian of the stencil w.r.t. the stencil's vertex positions.
    virtual MatrixMax12d
    compute_distance_hessian(const VectorMax12d& positions) const = 0;
};

} // namespace ipc