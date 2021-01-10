#pragma once

#include <ipc/collision_constraint.hpp>
#include <ipc/friction/relative_displacement.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

// C1 clamping
template <typename T>
inline T f0_SF(const T& x_squared, const double& epsv_times_h)
{
    double epsv_times_h_squared = epsv_times_h * epsv_times_h;
    if (x_squared >= epsv_times_h_squared) {
        return sqrt(x_squared);
    }
    return x_squared * (-sqrt(x_squared) / 3.0 + epsv_times_h)
        / (epsv_times_h_squared)
        + epsv_times_h / 3.0;
}

/// Derivative of f0_SF divided by the relative norm
template <typename T>
inline T f1_SF_div_relative_displacement_norm(
    const T& x_squared, const double& epsv_times_h)
{
    double epsv_times_h_squared = epsv_times_h * epsv_times_h;
    if (x_squared >= epsv_times_h_squared) {
        return 1 / sqrt(x_squared);
    }
    return (-sqrt(x_squared) + 2.0 * epsv_times_h) / epsv_times_h_squared;
}

template <typename T> inline T f2_SF(const T&, const double& epsv_times_h)
{
    // same for x_squared >= epsv_times_h * epsv_times_h for C1 clamped friction
    return T(-1 / (epsv_times_h * epsv_times_h));
}

///////////////////////////////////////////////////////////////////////////////

struct FrictionConstraint {
    /// @brief Barycentric coordinates of the closest point(s)
    Eigen::VectorX2d closest_point;

    /// @brief Tangent basis of the contact (max size 3×2)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 3, 2>
        tangent_basis;

    /// @brief Contact force magnitude
    double normal_force_magnitude;

    /// @brief Coefficient of friction
    double mu;

    virtual ~FrictionConstraint() {}

    virtual std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const = 0;

    virtual Eigen::VectorX12d compute_potential_gradient(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const = 0;

    virtual Eigen::MatrixXX12d compute_potential_hessian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h,
        const bool project_to_psd) const = 0;

protected:
    template <typename DerivedRelUi, typename T = typename DerivedRelUi::Scalar>
    T compute_potential_common(
        const Eigen::MatrixBase<DerivedRelUi>& rel_ui,
        const double epsv_times_h) const
    {
        return mu * normal_force_magnitude
            * f0_SF(
                   (rel_ui.transpose() * tangent_basis.cast<T>()).squaredNorm(),
                   epsv_times_h);
    }

    Eigen::VectorX2d compute_potential_gradient_common(
        const Eigen::VectorX3d& relative_displacement,
        double epsv_times_h) const;

    Eigen::MatrixXX12d compute_potential_hessian_common(
        const Eigen::VectorX3d& relative_displacement,
        const Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::ColMajor, 2, 12>&
            TT,
        const double epsv_times_h,
        bool project_to_psd,
        const int multiplicity = 1) const;
};

///////////////////////////////////////////////////////////////////////////////

struct VertexVertexFrictionConstraint : VertexVertexCandidate,
                                        FrictionConstraint {
    VertexVertexFrictionConstraint(long vertex0_index, long vertex1_index);
    VertexVertexFrictionConstraint(const VertexVertexCandidate& candidate);
    VertexVertexFrictionConstraint(const VertexVertexConstraint& constraint);

    std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex0_index, vertex1_index } };
    }

    template <typename T>
    Eigen::VectorX3<T> relative_displacement(const Eigen::MatrixX<T>& U) const
    {
        return point_point_relative_displacement(
            U.row(vertex0_index), U.row(vertex1_index));
    }

    template <typename T>
    T compute_potential(
        const Eigen::MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        Eigen::VectorX3<T> rel_u = relative_displacement(U);
        return multiplicity * compute_potential_common(rel_u, epsv_times_h);
    }

    Eigen::VectorX12d compute_potential_gradient(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const override;

    Eigen::MatrixXX12d compute_potential_hessian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h,
        const bool project_to_psd) const override;

    long multiplicity = 1;
};

///////////////////////////////////////////////////////////////////////////////

struct EdgeVertexFrictionConstraint : EdgeVertexCandidate, FrictionConstraint {
    EdgeVertexFrictionConstraint(long edge_index, long vertex_index);
    EdgeVertexFrictionConstraint(const EdgeVertexCandidate& constraint);
    EdgeVertexFrictionConstraint(const EdgeVertexConstraint& constraint);

    std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, E(edge_index, 0), E(edge_index, 1) } };
    }

    template <typename T>
    Eigen::VectorX3<T> relative_displacement(
        const Eigen::MatrixX<T>& U, const Eigen::MatrixXi& E) const
    {
        return point_edge_relative_displacement(
            U.row(vertex_index), //
            U.row(E(edge_index, 0)), U.row(E(edge_index, 1)),
            T(closest_point[0]));
    }

    template <typename T>
    T compute_potential(
        const Eigen::MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        Eigen::VectorX3<T> rel_u = relative_displacement(U, E);
        return multiplicity * compute_potential_common(rel_u, epsv_times_h);
    }

    Eigen::VectorX12d compute_potential_gradient(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const override;

    Eigen::MatrixXX12d compute_potential_hessian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h,
        const bool project_to_psd) const override;

    long multiplicity = 1;
};

///////////////////////////////////////////////////////////////////////////////

struct EdgeEdgeFrictionConstraint : EdgeEdgeCandidate, FrictionConstraint {
    EdgeEdgeFrictionConstraint(long edge0_index, long edge1_index);
    EdgeEdgeFrictionConstraint(const EdgeEdgeCandidate& constraint);
    EdgeEdgeFrictionConstraint(const EdgeEdgeConstraint& constraint);

    std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { E(edge0_index, 0), E(edge0_index, 1), //
                   E(edge1_index, 0), E(edge1_index, 1) } };
    }

    template <typename T>
    Eigen::VectorX3<T> relative_displacement(
        const Eigen::MatrixX<T>& U, const Eigen::MatrixXi& E) const
    {
        return edge_edge_relative_displacement(
            U.row(E(edge0_index, 0)), U.row(E(edge0_index, 1)),
            U.row(E(edge1_index, 0)), U.row(E(edge1_index, 1)),
            closest_point.cast<T>());
    }

    template <typename T>
    T compute_potential(
        const Eigen::MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        return compute_potential_common(
            relative_displacement(U, E), epsv_times_h);
    }

    Eigen::VectorX12d compute_potential_gradient(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const override;

    Eigen::MatrixXX12d compute_potential_hessian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h,
        const bool project_to_psd) const override;
};

///////////////////////////////////////////////////////////////////////////////

struct FaceVertexFrictionConstraint : FaceVertexCandidate, FrictionConstraint {
    FaceVertexFrictionConstraint(long face_index, long vertex_index);
    FaceVertexFrictionConstraint(const FaceVertexCandidate& constraint);
    FaceVertexFrictionConstraint(const FaceVertexConstraint& constraint);

    std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, //
                   F(face_index, 0), F(face_index, 1), F(face_index, 2) } };
    }

    template <typename T>
    Eigen::VectorX3<T> relative_displacement(
        const Eigen::MatrixX<T>& U, const Eigen::MatrixXi& F) const
    {
        return point_triangle_relative_displacement(
            U.row(vertex_index), U.row(F(face_index, 0)),
            U.row(F(face_index, 1)), U.row(F(face_index, 2)),
            closest_point.cast<T>());
    }

    template <typename T>
    T compute_potential(
        const Eigen::MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        return compute_potential_common(
            relative_displacement(U, F), epsv_times_h);
    }

    Eigen::VectorX12d compute_potential_gradient(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const override;

    Eigen::MatrixXX12d compute_potential_hessian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h,
        const bool project_to_psd) const override;
};

///////////////////////////////////////////////////////////////////////////////

struct FrictionConstraints {
    std::vector<VertexVertexFrictionConstraint> vv_constraints;
    std::vector<EdgeVertexFrictionConstraint> ev_constraints;
    std::vector<EdgeEdgeFrictionConstraint> ee_constraints;
    std::vector<FaceVertexFrictionConstraint> fv_constraints;

    size_t size() const;

    size_t num_constraints() const;

    void clear();

    FrictionConstraint& operator[](size_t idx);
    const FrictionConstraint& operator[](size_t idx) const;
};

} // namespace ipc
