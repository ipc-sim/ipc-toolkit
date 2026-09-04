#pragma once

#include "../adaptive_support.hpp"
#include "esp_primitives.hpp"
#include "vertex_matrix_view.hpp"

#include <ipc/barrier/barrier.hpp>
#include <ipc/esp/esp_parameters.hpp>
#include <ipc/math/math.hpp>
#include <ipc/utils/autodiff_types.hpp>

#include <array>

namespace ipc {

enum class ESPCollisionType : uint8_t {
    EDGE_VERTEX = 0,
    VERTEX_VERTEX = 1,
    FACE_VERTEX = 2,
    EDGE_EDGE = 3,
    EDGE_FACE = 4,
    FACE_FACE = 5
};

/// @brief Contact pair class for Geometric Contact Potential.
/// @note Unlike NormalCollision, ESPCollision has to be reconstructed whenever vertices change position
class ESPCollision {
public:
    static constexpr int MAX_VERT_3D = 20 * 2;
    static constexpr int ELEMENT_SIZE = 3 * MAX_VERT_3D;

    ESPCollision() = default;

    virtual ~ESPCollision() = default;

    /// @brief Name of the contact pair type
    virtual std::string name() const = 0;

    /// @brief Number of vertices involved times the dimension
    virtual int n_dofs() const = 0;

    /// @brief Contact pair type
    virtual ESPCollisionType type() const = 0;

    virtual std::array<index_t, 3> get_typed_hash() const = 0;

    /// @brief Get the number of vertices in the collision stencil.
    virtual int num_vertices() const = 0;

    /// @brief Get the number of vertices in primitive A's stencil.
    virtual size_t n_vertices_a() const = 0;

    /// @brief Get the number of vertices in primitive B's stencil.
    virtual size_t n_vertices_b() const = 0;

    /// @brief Get the vertex IDs of the collision stencil.
    /// @return The vertex IDs of the collision stencil. Size is always 4, but elements i > num_vertices() are -1.
    std::vector<index_t> vertex_ids() const;
    virtual index_t vertex_id(index_t i) const = 0;

    /// @brief Get the vertex attributes of the collision stencil.
    /// @param vertices Vertex attributes
    /// @return The vertex positions of the collision stencil. Size is always 4, but elements i > num_vertices() are NaN.
    Eigen::MatrixXd vertices(Eigen::ConstRef<Eigen::MatrixXd> vertices) const
    {
        const int DIM = vertices.cols();
        Eigen::MatrixXd stencil_vertices(num_vertices(), DIM);
        for (int i = 0; i < num_vertices(); i++) {
            stencil_vertices.row(i) = vertices.row(vertex_id(i));
        }

        return stencil_vertices;
    }

    /// @brief Select this stencil's DOF from the full matrix of DOF.
    /// @param X Full matrix of DOF (rowwise).
    /// @return This stencil's DOF.
    Eigen::VectorXd dof(Eigen::ConstRef<Eigen::MatrixXd> X) const;

    /// @brief Select this stencil's DOF from the full matrix of DOF.
    /// In 3D, some vertices may not be directly stored in the full matrix, e.g.
    /// face centers and edge-edge closest points.
    Eigen::VectorXd dof(VertexMatrixView<3> X_extended) const;

    /// @brief Select this stencil's DOF from the full 2D matrix of DOF (with a virtual vertex appended).
    Eigen::VectorXd dof(VertexMatrixView<2> X_extended) const;

    /// @brief Compute the distance of the stencil.
    /// @param vertices Collision mesh vertices
    /// @return Squared distance of the stencil.
    virtual double
    compute_distance(Eigen::ConstRef<Eigen::MatrixXd> vertices) const = 0;

    /// @brief Compute the value of the GCP potential
    virtual double operator()(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive = nullptr) const = 0;

    /// @brief Compute the gradient of the GCP potential wrt. vertices involved
    virtual VectorMax<double, ELEMENT_SIZE> gradient(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive = nullptr) const = 0;

    /// @brief Compute the Hessian of the GCP potential wrt. vertices involved
    virtual MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE> hessian(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive = nullptr) const = 0;

    bool operator==(const ESPCollision& other) const
    {
        return ((*this)[0] == other[0] && (*this)[1] == other[1]);
    }

    bool operator!=(const ESPCollision& other) const
    {
        return !(*this == other);
    }

    virtual index_t operator[](int idx) const = 0;

    virtual std::pair<index_t, index_t> get_hash() const = 0;

    virtual std::pair<double, double> operator_nearfar(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const NearFarBarrier* nf_barrier) const
    {
        return { 0.0, 0.0 };
    }

    virtual std::
        pair<VectorMax<double, ELEMENT_SIZE>, VectorMax<double, ELEMENT_SIZE>>
        gradient_nearfar(
            Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
            const ESPParameters& params,
            const AdaptiveSupport* adaptive,
            const NearFarBarrier* nf_barrier) const
    {
        VectorMax<double, ELEMENT_SIZE> zero =
            VectorMax<double, ELEMENT_SIZE>::Zero(positions.size());
        return { zero, zero };
    }

    virtual std::pair<
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>,
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>>
    hessian_nearfar(
        Eigen::ConstRef<VectorMax<double, ELEMENT_SIZE>> positions,
        const ESPParameters& params,
        const AdaptiveSupport* adaptive,
        const NearFarBarrier* nf_barrier) const
    {
        int n = positions.size();
        MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE> zero =
            MatrixMax<double, ELEMENT_SIZE, ELEMENT_SIZE>::Zero(n, n);
        return { zero, zero };
    }

    double weight = 1;
};

} // namespace ipc