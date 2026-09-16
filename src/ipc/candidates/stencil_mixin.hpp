#pragma once

#include <ipc/config.hpp>
#include <ipc/ccd/default_narrow_phase_ccd.hpp>
#include <ipc/geometry/normal.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>

#include <array>
#include <cassert>
#include <limits>
#include <ostream>

namespace ipc {

/// @brief Stencil math that does not depend on the stencil type.
///
/// Split out of StencilMixin so these are compiled once rather than once per
/// stencil type.
class StencilMath {
public:
    /// @brief The maximum number of vertices in a collision stencil.
    static constexpr int STENCIL_SIZE = 4;

    /// @brief Compute diag((∂t/∂x)(∂t/∂x)ᵀ) efficiently (Eq. 11 of the paper).
    ///
    /// Result is [c₀², c₀², c₀², c₁², c₁², c₁², ...] (each cᵢ² repeated dim
    /// times). This is used for computing diag(∂²b/∂x²) without the full
    /// Hessian.
    ///
    /// @param coeffs The coefficients cᵢ (from compute_coefficients).
    /// @param d The spatial dimension (2 or 3).
    /// @return The diagonal of (∂t/∂x)(∂t/∂x)ᵀ as a vector of size ndof.
    static VectorMax12d
    diag_distance_vector_outer(Eigen::ConstRef<VectorMax4d> coeffs, const int d)
    {
        const int n = coeffs.size();
        VectorMax12d diag(n * d);
        for (int i = 0; i < n; i++) {
            const double ci2 = coeffs[i] * coeffs[i];
            diag.segment(i * d, d).setConstant(ci2);
        }
        return diag;
    }

    /// @brief Compute diag((∂t/∂x · t)(∂t/∂x · t)ᵀ) efficiently (Eq. 12).
    ///
    /// Result is element-wise square of [c₀tᵀ, c₁tᵀ, ..., cₙtᵀ].
    /// This is used for computing diag(∂²b/∂x²) without the full Hessian.
    ///
    /// @param coeffs The coefficients cᵢ.
    /// @param distance_vector The distance vector t.
    /// @return The diagonal of (∂t/∂x·t)(∂t/∂x·t)ᵀ as a vector of size ndof.
    static VectorMax12d diag_distance_vector_t_outer(
        Eigen::ConstRef<VectorMax4d> coeffs,
        Eigen::ConstRef<VectorMax3d> distance_vector)
    {
        const int n = coeffs.size();
        const int d = distance_vector.size();
        const VectorMax3d t2 = distance_vector.array().square();
        VectorMax12d diag(n * d);
        for (int i = 0; i < n; i++) {
            diag.segment(i * d, d).array() = (coeffs[i] * coeffs[i]) * t2;
        }
        return diag;
    }

    /// @brief Compute pᵀ(∂t/∂x) efficiently as ∑ cᵢ pᵢ (Eqs. 13-14).
    ///
    /// Given p = [p₀, p₁, ..., pₙ]ᵀ where pᵢ ∈ ℝ^dim, this computes
    /// pᵀ(∂t/∂x) = ∑ cᵢ pᵢ which is a dim-dimensional vector.
    ///
    /// @param coeffs The coefficients cᵢ.
    /// @param p A vector of size ndof (the direction for the quadratic form).
    /// @param d The spatial dimension (2 or 3).
    /// @return pᵀ(∂t/∂x) as a dim-dimensional vector.
    static VectorMax3d contract_distance_vector_jacobian(
        Eigen::ConstRef<VectorMax4d> coeffs,
        Eigen::ConstRef<VectorMax12d> p,
        const int d)
    {
        const int n = coeffs.size();
        VectorMax3d result = VectorMax3d::Zero(d);
        for (int i = 0; i < n; i++) {
            result += coeffs[i] * p.segment(d * i, d);
        }
        return result;
    }
};

/// @brief The stencil operations that are written once in terms of a few
/// primitives every stencil provides.
///
/// A stencil supplies `num_vertices()`, `vertex_ids()`, the four
/// `positions`-based quantities (distance, its gradient and Hessian, and the
/// coefficients), `ccd()`, and the two unnormalized-normal quantities. This
/// mixin builds everything else on top of them.
///
/// It is deliberately empty: no data members and no virtual functions. A
/// stencil that derives from it therefore pays nothing in size (empty base
/// optimization), and stays trivially copyable if its own members are, which
/// is what lets the broad phase move millions of candidates with a memcpy.
/// `CollisionStencil` derives from it too, and there the calls below land on
/// its virtual functions.
///
/// @tparam Derived The stencil type deriving from this.
template <class Derived> class StencilMixin : public StencilMath {
protected:
    constexpr const Derived& self() const
    {
        return static_cast<const Derived&>(*this);
    }

public:
    /// @brief Get the dimension of the collision stencil.
    /// @param ndof Number of degrees of freedom in the stencil.
    /// @return The dimension of the collision stencil.
    int dim(const int ndof) const
    {
        assert(ndof % self().num_vertices() == 0);
        return ndof / self().num_vertices();
    }

    /// @brief Get the vertex attributes of the collision stencil.
    /// @param vertices Vertex attributes
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @return The vertex positions of the collision stencil. Elements i > num_vertices() are NaN.
    std::array<VectorMax3d, STENCIL_SIZE> vertices(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        constexpr double NaN = std::numeric_limits<double>::signaling_NaN();

        const auto vertex_ids = self().vertex_ids(edges, faces);

        std::array<VectorMax3d, STENCIL_SIZE> stencil_vertices;
        for (int i = 0; i < STENCIL_SIZE; i++) {
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
        const int n = self().num_vertices();
        VectorMax12d x(n * dim);
        const auto idx = self().vertex_ids(edges, faces);
        for (int i = 0; i < n; i++) {
            x.segment(i * dim, dim) = X.row(idx[i]);
        }
        return x;
    }

    // -- Mesh-based overloads of the positions-based quantities --------------

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
        return self().compute_distance(dof(vertices, edges, faces));
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
        return self().compute_distance_gradient(dof(vertices, edges, faces));
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
        return self().compute_distance_hessian(dof(vertices, edges, faces));
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
        return self().compute_coefficients(dof(vertices, edges, faces));
    }

    /// @brief Compute the distance vector using the mesh vertices.
    /// @param vertices Collision mesh vertices.
    /// @param edges Collision mesh edges.
    /// @param faces Collision mesh faces.
    /// @return The distance vector (dim-dimensional).
    VectorMax3d compute_distance_vector(
        Eigen::ConstRef<Eigen::MatrixXd> vertices,
        Eigen::ConstRef<Eigen::MatrixXi> edges,
        Eigen::ConstRef<Eigen::MatrixXi> faces) const
    {
        return compute_distance_vector(dof(vertices, edges, faces));
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

    /// @brief Compute the distance vector of the stencil: t = ∑ cᵢ xᵢ.
    ///
    /// The distance vector is the vector between the closest points on the
    /// collision primitives. Its squared norm equals the squared distance.
    ///
    /// @param positions Stencil's vertex positions.
    /// @return The distance vector (dim-dimensional, i.e., 2D or 3D).
    VectorMax3d
    compute_distance_vector(Eigen::ConstRef<VectorMax12d> positions) const
    {
        VectorMax4d _; // Unused output
        return compute_distance_vector(positions, _);
    }

    /// @brief Compute the distance vector and the coefficients together.
    /// @param positions Stencil's vertex positions.
    /// @param[out] coeffs The computed coefficients cᵢ.
    /// @return The distance vector.
    VectorMax3d compute_distance_vector(
        Eigen::ConstRef<VectorMax12d> positions, VectorMax4d& coeffs) const
    {
        const int n = self().num_vertices();
        const int d = dim(positions.size());
        coeffs = self().compute_coefficients(positions);
        VectorMax3d dv = VectorMax3d::Zero(d);
        for (int i = 0; i < n; i++) {
            dv += coeffs[i] * positions.segment(d * i, d);
        }
        return dv;
    }

    /// @brief Compute the Jacobian of the distance vector w.r.t. positions: ∂t/∂x = [c₀I c₁I ... cₙI]ᵀ.
    ///
    /// Since ∂t/∂x has a very simple structure (block-diagonal with scalar
    /// coefficients times identity), many operations can be done without
    /// forming this matrix. This method is provided for completeness and
    /// verification.
    ///
    /// @param positions Stencil's vertex positions.
    /// @return The Jacobian ∂t/∂x ∈ ℝ^{ndof × dim}.
    MatrixMax<double, 12, 3> compute_distance_vector_jacobian(
        Eigen::ConstRef<VectorMax12d> positions) const
    {
        const int n = self().num_vertices();
        const int d = dim(positions.size());
        const int ndof = n * d;
        const VectorMax4d c = self().compute_coefficients(positions);
        // ∂t/∂x is ndof × dim
        // Each block row i is cᵢ * I_{d×d}
        MatrixMax<double, 12, 3> J;
        J.setZero(ndof, d);
        for (int i = 0; i < n; i++) {
            J.block(i * d, 0, d, d).diagonal().array() = c[i];
        }
        return J;
    }

    /// @brief Compute the normal of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @param flip_if_negative If true, flip the normal if the point is on the negative side.
    /// @param sign If not nullptr, set to the sign of the normal before any flipping.
    /// @return Normal of the stencil.
    VectorMax3d compute_normal(
        Eigen::ConstRef<VectorMax12d> positions,
        bool flip_if_negative = true,
        double* sign = nullptr) const
    {
        const int dim = this->dim(positions.size());

        VectorMax3d n =
            self().compute_unnormalized_normal(positions).normalized();

        if (sign != nullptr) {
            *sign =
                (positions.head(dim) - positions.tail(dim)).dot(n) < 0 ? -1 : 1;
        }

        // Flip the normal if the point is on the negative side.
        // Any point on the second object will do, so we use the last point.
        if (flip_if_negative
            && (positions.head(dim) - positions.tail(dim)).dot(n) < 0) {
            n *= -1;
        }

        return n;
    }

    /// @brief Compute the Jacobian of the normal of the stencil.
    /// @param positions Stencil's vertex positions.
    /// @param flip_if_negative If true, flip the normal if the point is on the negative side.
    /// @return Jacobian of the normal of the stencil.
    MatrixMax<double, 3, 12> compute_normal_jacobian(
        Eigen::ConstRef<VectorMax12d> positions,
        bool flip_if_negative = true) const
    {
        const int dim = this->dim(positions.size());

        const VectorMax3d n = self().compute_unnormalized_normal(positions);

        MatrixMax<double, 3, 12> dn = normalization_jacobian(n)
            * self().compute_unnormalized_normal_jacobian(positions);

        if (flip_if_negative
            && (positions.head(dim) - positions.tail(dim)).dot(n) < 0) {
            dn *= -1;
        }

        return dn;
    }

    /// @brief Write the CCD query to a stream.
    /// @param out Stream to write to.
    /// @param vertices_t0 Stencil vertices at the start of the time step.
    /// @param vertices_t1 Stencil vertices at the end of the time step.
    /// @return The stream.
    std::ostream& write_ccd_query(
        std::ostream& out,
        Eigen::ConstRef<VectorMax12d> vertices_t0,
        Eigen::ConstRef<VectorMax12d> vertices_t1) const
    {
        assert(vertices_t0.size() == vertices_t1.size());

        const int n = self().num_vertices();
        const int dim = vertices_t0.size() / n;
        assert(vertices_t0.size() % n == 0);

        for (int i = 0; i < n; i++) {
            out << vertices_t0.segment(dim * i, dim)
                       .transpose()
                       .format(OBJ_VERTEX_FORMAT);
        }

        for (int i = 0; i < n; i++) {
            out << vertices_t1.segment(dim * i, dim)
                       .transpose()
                       .format(OBJ_VERTEX_FORMAT);
        }

        return out;
    }
};

} // namespace ipc
