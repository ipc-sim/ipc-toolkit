#pragma once

#include "primitive.hpp"

#include <ipc/smooth_contact/distance/mollifier.hpp>

namespace ipc {
class Edge3 : public Primitive {
public:
    static constexpr int N_CORE_POINTS = 2;
    static constexpr int DIM = 3;
    static constexpr int MAX_SIZE = N_EDGE_NEIGHBORS_3D * DIM;

    /// @brief Construct a smooth edge3 primitive.
    /// @param id Primitive ID
    /// @param mesh Collision mesh containing the edge
    /// @param vertices Vertex positions of the mesh
    /// @param d Vector from closest point on the edge to the point outside of the edge
    /// @param params Smooth contact parameters
    Edge3(
        const index_t id,
        const CollisionMesh& mesh,
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<VectorMax3d> d,
        const SmoothContactParameters& params);

    /// @brief Get the number of vertices (edge endpoints + face-opposite vertices)
    int n_vertices() const override { return m_vertex_ids.size(); }

    /// @brief Get the number of DOFs (dim per vertex)
    int n_dofs() const override { return n_vertices() * DIM; }

    /// @brief Number of adjacent face neighbors (opposite vertices)
    int n_face_neighbors() const { return n_neighbors; }

    // assume the following functions are only called if active

    /// @brief Potential value of the smooth edge3 term.
    /// The order of DOFs in x is [e0(3), e1(3), f0(3), f1(3), ...] where fi are
    /// the opposite vertices of the face neighbors.
    /// @param d The closest direction from the edge to the point outside of the edge
    /// @param x The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1(3), ...]
    /// @return The potential value of the smooth edge3 term
    double potential(
        Eigen::ConstRef<Eigen::Vector3d> d,
        Eigen::ConstRef<VectorMax<double, MAX_SIZE>> x) const;

    // derivatives including wrt. d (the closest direction) in front

    /// @brief Gradient with respect to the edge vertices and face-opposite vertices,
    /// as well as the closest direction d. The order of DOFs in the output is
    /// [d(3), e0(3), e1(3), f0(3), f1(3), ...] where fi are the opposite
    /// vertices of the face neighbors.
    /// @param d The closest direction from the edge to the point outside of the edge
    /// @param x The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1(3), ...]
    /// @return The gradient of the smooth edge3 term with respect to [d, e0, e1, f0, f1, ...]
    VectorMax<double, MAX_SIZE + DIM> grad(
        Eigen::ConstRef<Eigen::Vector3d> d,
        Eigen::ConstRef<VectorMax<double, MAX_SIZE>> x) const;

    /// @brief Hessian with respect to the edge vertices and face-opposite vertices,
    /// as well as the closest direction d. The order of DOFs in the output is
    /// [d(3), e0(3), e1(3), f0(3), f1(3), ...] where fi are the opposite
    /// vertices of the face neighbors
    /// @param d The closest direction from the edge to the point outside of the edge
    /// @param x The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1(3), ...]
    /// @return The Hessian of the smooth edge3 term with respect to [d, e0, e1, f0, f1, ...]
    MatrixMax<double, MAX_SIZE + DIM, MAX_SIZE + DIM> hessian(
        Eigen::ConstRef<Eigen::Vector3d> d,
        Eigen::ConstRef<VectorMax<double, MAX_SIZE>> x) const;

    /// @brief Templated implementation for autodiff verification
    /// @tparam T Scalar type for autodiff (e.g., ADGrad or ADHessian)
    /// @tparam n_verts Number of vertices (edge endpoints + face-opposite vertices).
    /// @param X The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1(3), ...]
    /// @param direction The closest direction from the edge to the point outside of the edge
    /// @return The potential value of the smooth edge3 term
    template <typename T, int n_verts = Eigen::Dynamic>
    T smooth_edge3_term(
        Eigen::ConstRef<Eigen::Matrix<T, n_verts, 3>> X,
        Eigen::ConstRef<Eigen::RowVector3<T>> direction) const;

    /// @brief Templated gradient for autodiff verification
    /// @tparam T Scalar type for autodiff (e.g., ADGrad or ADHessian)
    /// @param direction The closest direction from the edge to the point outside of the edge
    /// @param X The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1(3), ...]
    /// @return The gradient of the smooth edge3 term with respect to [direction, e0, e1, f0, f1, ...]
    GradientType<Eigen::Dynamic> smooth_edge3_term_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> direction,
        Eigen::ConstRef<Eigen::MatrixX3d> X) const;

    /// @brief Templated Hessian for autodiff verification
    /// @tparam T Scalar type for autodiff (e.g., ADGrad or ADHessian)
    /// @param direction The closest direction from the edge to the point outside of the edge
    /// @param X The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1, ...]
    /// @return The Hessian of the smooth edge3 term with respect to [direction, e0, e1, f0, f1, ...]
    HessianType<Eigen::Dynamic> smooth_edge3_term_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> direction,
        Eigen::ConstRef<Eigen::MatrixX3d> X) const;

    /// @brief Gradient of the normal term with respect to the edge vertices and face-opposite vertices, as well as the closest direction d.
    /// @param direction The normalized closest direction from the point outside of the edge to the edge
    ///        (i.e., direction = -d.normalized(), where d is the vector from
    ///        the edge to the outside point)
    /// @param X The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1(3), ...]
    /// @param alpha The alpha parameter for the normal term heaviside function
    /// @param beta The beta parameter for the normal term heaviside function
    /// @return The gradient of the normal term with respect to [d, e0, e1, f0, f1, ...], where d is defined as above
    GradientType<Eigen::Dynamic> smooth_edge3_normal_term_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> direction,
        Eigen::ConstRef<Eigen::MatrixX3d> X,
        const double alpha,
        const double beta) const;

    /// @brief Hessian of the normal term with respect to the edge vertices and face-opposite vertices, as well as the closest direction d.
    /// @param direction The normalized closest direction from the point outside of the edge to the edge
    ///        (i.e., direction = -d.normalized(), where d is the vector from
    ///        the edge to the outside point)
    /// @param X The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1(3), ...]
    /// @param alpha The alpha parameter for the normal term heaviside function
    /// @param beta The beta parameter for the normal term heaviside function
    /// @return The Hessian of the normal term with respect to [d, e0, e1, f0, f1, ...], where d is defined as above
    HessianType<Eigen::Dynamic> smooth_edge3_normal_term_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> direction,
        Eigen::ConstRef<Eigen::MatrixX3d> X,
        const double alpha,
        const double beta) const;

    /// @brief Gradient of the tangent term with respect to the edge vertices and face-opposite vertices, as well as the closest direction d.
    /// The order of DOFs in the output is [d(3), e0(3), e1(3), f0(3), f1(3),
    /// ...] where fi are the opposite vertices of the face neighbors.
    /// @param dn The normalized, already-negated closest direction (from the point outside of the edge to the edge), i.e., -d / ||d||, matching the -dn·t/|t| convention.
    /// @param tangents The tangent directions for each face neighbor, computed as point_line_closest_point_direction(fi, e0, e1) for each face neighbor. The order of rows is the same as the order of face neighbors in faces.
    /// @param alpha The alpha parameter for the tangent term heaviside function
    /// @param beta The beta parameter for the tangent term heaviside function
    /// @return The gradient of the tangent term with respect to [d, e0, e1, f0, f1, ...]
    GradientType<Eigen::Dynamic> smooth_edge3_tangent_term_gradient(
        Eigen::ConstRef<Eigen::RowVector3d> dn,
        Eigen::ConstRef<Eigen::MatrixX3d> tangents,
        const double alpha,
        const double beta) const;

    /// @brief Hessian of the tangent term with respect to the edge vertices and face-opposite vertices, as well as the closest direction d.
    /// @param dn The normalized, already-negated closest direction (from the point outside of the edge to the edge), i.e., -d / ||d||, matching the -dn·t/|t| convention.
    /// @param tangents The tangent directions for each face neighbor, computed as point_line_closest_point_direction(fi, e0, e1) for each face neighbor. The order of rows is the same as the order of face neighbors in faces.
    /// @param alpha The alpha parameter for the tangent term heaviside function
    /// @param beta The beta parameter for the tangent term heaviside function
    /// @return The Hessian of the tangent term with respect to [d, e0, e1, f0, f1, ...]
    HessianType<Eigen::Dynamic> smooth_edge3_tangent_term_hessian(
        Eigen::ConstRef<Eigen::RowVector3d> dn,
        Eigen::ConstRef<Eigen::MatrixX3d> tangents,
        const double alpha,
        const double beta) const;

private:
    /// @brief Check if the smooth edge3 term is active (i.e., if the tangent and normal terms are not trivially 1)
    /// @param X The positions of the edge vertices and face-opposite vertices, in the order [e0(3), e1(3), f0(3), f1(3), ...]
    /// @param direction The closest direction from the point outside of the edge to the edge
    ///        (i.e., direction = -d.normalized(), where d is the vector from
    ///        the edge to the outside point)
    /// @return True if the smooth edge3 term is active, false if it is trivially 1
    bool smooth_edge3_term_type(
        Eigen::ConstRef<Eigen::MatrixX3d> X,
        Eigen::ConstRef<Eigen::RowVector3d> direction);

    /// @brief Number of face neighbors (opposite vertices)
    int n_neighbors;

    /// @brief Orientation types for tangent and normal terms, indexed by face neighbor
    OrientationTypes otypes;

    /// @brief Local face connectivity: each row is [local_e0, local_e1, local_fi]
    /// where local_e0/e1 are the edge vertex local indices (0 and 1) and
    /// local_fi is the local index of the opposite vertex.
    /// Size: n_neighbors × 3
    Eigen::MatrixX3i faces;

    /// @brief Whether the edge is orientable. If false, the normal term is not applied.
    bool orientable;
    double m_rest_sq_length; ///< Rest-shape squared edge length (quadrature
                             ///< weight, constant)
};

} // namespace ipc