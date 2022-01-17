#pragma once

#include <ipc/collision_constraint.hpp>
#include <ipc/friction/relative_displacement.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct FrictionConstraint {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /// @brief Barycentric coordinates of the closest point(s)
    VectorMax2d closest_point;

    /// @brief Tangent basis of the contact (max size 3×2)
    MatrixMax<double, 3, 2> tangent_basis;

    /// @brief Contact force magnitude
    double normal_force_magnitude;

    /// @brief Coefficient of friction
    double mu;

    virtual ~FrictionConstraint() { }

    virtual std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const = 0;

    virtual VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const;

    virtual MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h,
        const bool project_hessian_to_psd) const;

    virtual double compute_normal_force_magnitude(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const = 0;

    virtual VectorMax12d compute_normal_force_magnitude_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const = 0;

    virtual MatrixMax<double, 3, 2> compute_tangent_basis(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual MatrixMax<double, 3, 24> compute_tangent_basis_jacobian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual VectorMax3d relative_displacement(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual MatrixMax<double, 3, 12> relative_displacement_jacobian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

protected:
    virtual int get_multiplicity() const { return 1; };

    template <typename DerivedRelUi, typename T = typename DerivedRelUi::Scalar>
    T compute_potential_common(
        const Eigen::MatrixBase<DerivedRelUi>& rel_ui,
        const double epsv_times_h) const
    {
        // u is the relative displacement in the tangential space
        const VectorMax2d u = tangent_basis.transpose().cast<T>() * rel_ui;
        return mu * normal_force_magnitude * f0_SF(u.norm(), epsv_times_h);
    }
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
    T compute_potential(
        const MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        VectorMax3<T> rel_u = relative_displacement_T(U);
        return multiplicity * compute_potential_common(rel_u, epsv_times_h);
    }

    double compute_normal_force_magnitude(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const override;

    VectorMax12d compute_normal_force_magnitude_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const override;

    MatrixMax<double, 3, 2> compute_tangent_basis(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax<double, 3, 24> compute_tangent_basis_jacobian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax3d relative_displacement(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return relative_displacement_T(U);
    }

    MatrixMax<double, 3, 12> relative_displacement_jacobian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    long multiplicity = 1;

protected:
    int get_multiplicity() const override { return multiplicity; };

private:
    template <typename T>
    VectorMax3<T> relative_displacement_T(const MatrixX<T>& U) const
    {
        return point_point_relative_displacement(
            U.row(vertex0_index), U.row(vertex1_index));
    }
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
    T compute_potential(
        const MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        VectorMax3<T> rel_u = relative_displacement_T(U, E);
        return multiplicity * compute_potential_common(rel_u, epsv_times_h);
    }

    double compute_normal_force_magnitude(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const override;

    VectorMax12d compute_normal_force_magnitude_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const override;

    MatrixMax<double, 3, 2> compute_tangent_basis(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax<double, 3, 24> compute_tangent_basis_jacobian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax3d relative_displacement(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return relative_displacement_T(U, E);
    }

    MatrixMax<double, 3, 12> relative_displacement_jacobian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    long multiplicity = 1;

protected:
    int get_multiplicity() const override { return multiplicity; };

private:
    template <typename T>
    VectorMax3<T>
    relative_displacement_T(const MatrixX<T>& U, const Eigen::MatrixXi& E) const
    {
        return point_edge_relative_displacement(
            U.row(vertex_index), //
            U.row(E(edge_index, 0)), U.row(E(edge_index, 1)),
            T(closest_point[0]));
    }
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
    T compute_potential(
        const MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        return compute_potential_common(
            relative_displacement_T(U, E), epsv_times_h);
    }

    double compute_normal_force_magnitude(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const override;

    VectorMax12d compute_normal_force_magnitude_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const override;

    MatrixMax<double, 3, 2> compute_tangent_basis(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax<double, 3, 24> compute_tangent_basis_jacobian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax3d relative_displacement(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return relative_displacement_T(U, E);
    }

    MatrixMax<double, 3, 12> relative_displacement_jacobian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

private:
    template <typename T>
    VectorMax3<T>
    relative_displacement_T(const MatrixX<T>& U, const Eigen::MatrixXi& E) const
    {
        return edge_edge_relative_displacement(
            U.row(E(edge0_index, 0)), U.row(E(edge0_index, 1)),
            U.row(E(edge1_index, 0)), U.row(E(edge1_index, 1)),
            closest_point.cast<T>());
    }
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
    T compute_potential(
        const MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        return compute_potential_common(
            relative_displacement_T(U, F), epsv_times_h);
    }

    double compute_normal_force_magnitude(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const override;

    VectorMax12d compute_normal_force_magnitude_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const override;

    MatrixMax<double, 3, 2> compute_tangent_basis(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax<double, 3, 24> compute_tangent_basis_jacobian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax3d relative_displacement(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override
    {
        return relative_displacement_T(U, F);
    }

    MatrixMax<double, 3, 12> relative_displacement_jacobian(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

private:
    template <typename T>
    VectorMax3<T>
    relative_displacement_T(const MatrixX<T>& U, const Eigen::MatrixXi& F) const
    {
        return point_triangle_relative_displacement(
            U.row(vertex_index), U.row(F(face_index, 0)),
            U.row(F(face_index, 1)), U.row(F(face_index, 2)),
            closest_point.cast<T>());
    }
};

///////////////////////////////////////////////////////////////////////////////

struct FrictionConstraints {
    template <typename T>
    using aligned_vector = std::vector<T, Eigen::aligned_allocator<T>>;

    aligned_vector<VertexVertexFrictionConstraint> vv_constraints;
    aligned_vector<EdgeVertexFrictionConstraint> ev_constraints;
    aligned_vector<EdgeEdgeFrictionConstraint> ee_constraints;
    aligned_vector<FaceVertexFrictionConstraint> fv_constraints;

    size_t size() const;

    size_t num_constraints() const;

    bool empty() const;

    void clear();

    FrictionConstraint& operator[](size_t idx);
    const FrictionConstraint& operator[](size_t idx) const;
};

} // namespace ipc
