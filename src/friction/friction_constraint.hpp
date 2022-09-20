#pragma once

#include <ipc/collisions/constraints.hpp>
#include <ipc/friction/relative_displacement.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <ipc/config.hpp>

namespace ipc {

struct FrictionConstraint {
    /// @brief Contact force magnitude
    double normal_force_magnitude;

    /// @brief Coefficient of friction
    double mu;

    double weight = 1;
    /// @brief Gradient of weight with respect to all DOF
    Eigen::SparseVector<double> weight_gradient;

    /// @brief Barycentric coordinates of the closest point(s)
    VectorMax2d closest_point;

    /// @brief Tangent basis of the contact (max size 3×2)
    MatrixMax<double, 3, 2> tangent_basis;

    virtual ~FrictionConstraint() { }

    virtual int num_vertices() const = 0;
    virtual std::array<long, 4> vertex_indices(
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

    virtual VectorMax12d compute_force(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) const //< whether to not multiply by mu
    {
        return compute_force(
            X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U, E, F, dhat,
            barrier_stiffness, epsv_times_h, dmin, no_mu);
    }

    virtual VectorMax12d compute_force(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) const; //< whether to not multiply by mu

    enum class DiffWRT { X, Ut, U };

    virtual MatrixMax12d compute_force_jacobian(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const DiffWRT wrt,
        const double dmin = 0) const
    {
        return compute_force_jacobian(
            X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U, E, F, dhat,
            barrier_stiffness, epsv_times_h, wrt, dmin);
    }

    virtual MatrixMax12d compute_force_jacobian(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const DiffWRT wrt,
        const double dmin = 0) const;

protected:
    int dim() const { return tangent_basis.rows(); }
    int ndof() const { return dim() * num_vertices(); };

    void init(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness,
        const double dmin);

    template <typename T>
    VectorMax12<T> select_dofs(
        const MatrixX<T>& X,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const
    {
        VectorMax12<T> x(ndof());
        const std::array<long, 4> idx = vertex_indices(E, F);
        for (int i = 0; i < num_vertices(); i++) {
            x.segment(i * dim(), dim()) = X.row(idx[i]);
        }
        return x;
    }

    virtual double compute_distance(const VectorMax12d& x) const = 0;
    virtual VectorMax12d
    compute_distance_gradient(const VectorMax12d& x) const = 0;

    virtual double compute_normal_force_magnitude(
        const VectorMax12d& x,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const;

    virtual VectorMax12d compute_normal_force_magnitude_gradient(
        const VectorMax12d& x,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const;

    virtual MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& x) const = 0;

    virtual MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& x) const = 0;

    virtual VectorMax2d compute_closest_point(const VectorMax12d& x) const = 0;

    virtual MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& x) const = 0;

    virtual VectorMax3d relative_displacement(const VectorMax12d& u) const = 0;

    virtual MatrixMax<double, 3, 12> relative_displacement_matrix() const
    {
        return relative_displacement_matrix(closest_point);
    }

    virtual MatrixMax<double, 3, 12>
    relative_displacement_matrix(const VectorMax2d& closest_point) const = 0;

    virtual MatrixMax<double, 6, 12> relative_displacement_matrix_jacobian(
        const VectorMax2d& closest_point) const = 0;

    template <typename DerivedRelUi, typename T = typename DerivedRelUi::Scalar>
    T compute_potential_common(
        const Eigen::MatrixBase<DerivedRelUi>& rel_ui,
        const double epsv_times_h) const
    {
        // u is the relative displacement in the tangential space
        const VectorMax2d u = tangent_basis.transpose().cast<T>() * rel_ui;
        return weight * mu * normal_force_magnitude
            * f0_SF(u.norm(), epsv_times_h);
    }
};

///////////////////////////////////////////////////////////////////////////////

struct VertexVertexFrictionConstraint : VertexVertexCandidate,
                                        FrictionConstraint {
    VertexVertexFrictionConstraint(long vertex0_index, long vertex1_index);
    VertexVertexFrictionConstraint(const VertexVertexCandidate& candidate);
    VertexVertexFrictionConstraint(const VertexVertexConstraint& constraint);
    VertexVertexFrictionConstraint(
        const VertexVertexConstraint& constraint,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness)
        : VertexVertexFrictionConstraint(constraint)
    {
        FrictionConstraint::init(
            V, E, F, dhat, barrier_stiffness, constraint.minimum_distance);
    }

    int num_vertices() const override { return 2; }
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex0_index, vertex1_index, -1, -1 } };
    }

    template <typename T>
    T compute_potential(
        const MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        return compute_potential_common(
            relative_displacement_T(select_dofs(U, E, F)), epsv_times_h);
    }

protected:
    virtual double compute_distance(const VectorMax12d& x) const override;
    virtual VectorMax12d
    compute_distance_gradient(const VectorMax12d& x) const override;

    MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& x) const override;

    MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& x) const override;

    VectorMax2d compute_closest_point(const VectorMax12d& x) const override;

    MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& x) const override;

    VectorMax3d relative_displacement(const VectorMax12d& u) const override
    {
        return relative_displacement_T(u);
    }

    using FrictionConstraint::relative_displacement_matrix;

    MatrixMax<double, 3, 12> relative_displacement_matrix(
        const VectorMax2d& closest_point) const override;

    MatrixMax<double, 6, 12> relative_displacement_matrix_jacobian(
        const VectorMax2d& closest_point) const override;

private:
    template <typename T>
    VectorMax3<T> relative_displacement_T(const VectorMax12<T>& u) const
    {
        assert(u.size() == ndof());
        return point_point_relative_displacement(u.head(dim()), u.tail(dim()));
    }
};

///////////////////////////////////////////////////////////////////////////////

struct EdgeVertexFrictionConstraint : EdgeVertexCandidate, FrictionConstraint {
    EdgeVertexFrictionConstraint(long edge_index, long vertex_index);
    EdgeVertexFrictionConstraint(const EdgeVertexCandidate& constraint);
    EdgeVertexFrictionConstraint(const EdgeVertexConstraint& constraint);
    EdgeVertexFrictionConstraint(
        const EdgeVertexConstraint& constraint,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness)
        : EdgeVertexFrictionConstraint(constraint)
    {
        FrictionConstraint::init(
            V, E, F, dhat, barrier_stiffness, constraint.minimum_distance);
    }

    int num_vertices() const override { return 3; }
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, E(edge_index, 0), E(edge_index, 1), -1 } };
    }

    template <typename T>
    T compute_potential(
        const MatrixX<T>& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const
    {
        return compute_potential_common(
            relative_displacement_T(select_dofs(U, E, F)), epsv_times_h);
    }

protected:
    virtual double compute_distance(const VectorMax12d& x) const override;
    virtual VectorMax12d
    compute_distance_gradient(const VectorMax12d& x) const override;

    MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& x) const override;

    MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& x) const override;

    VectorMax2d compute_closest_point(const VectorMax12d& x) const override;

    MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& x) const override;

    VectorMax3d relative_displacement(const VectorMax12d& u) const override
    {
        return relative_displacement_T(u);
    }

    using FrictionConstraint::relative_displacement_matrix;

    MatrixMax<double, 3, 12> relative_displacement_matrix(
        const VectorMax2d& closest_point) const override;

    MatrixMax<double, 6, 12> relative_displacement_matrix_jacobian(
        const VectorMax2d& closest_point) const override;

private:
    template <typename T>
    VectorMax3<T> relative_displacement_T(const VectorMax12<T>& u) const
    {
        assert(u.size() == ndof());
        return point_edge_relative_displacement(
            u.head(dim()), u.segment(dim(), dim()), u.tail(dim()),
            T(closest_point[0]));
    }
};

///////////////////////////////////////////////////////////////////////////////

struct EdgeEdgeFrictionConstraint : EdgeEdgeCandidate, FrictionConstraint {
    EdgeEdgeFrictionConstraint(long edge0_index, long edge1_index);
    EdgeEdgeFrictionConstraint(const EdgeEdgeCandidate& constraint);
    EdgeEdgeFrictionConstraint(const EdgeEdgeConstraint& constraint);
    EdgeEdgeFrictionConstraint(
        const EdgeEdgeConstraint& constraint,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness)
        : EdgeEdgeFrictionConstraint(constraint)
    {
        FrictionConstraint::init(
            V, E, F, dhat, barrier_stiffness, constraint.minimum_distance);
    }

    int num_vertices() const override { return 4; }
    std::array<long, 4> vertex_indices(
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
            relative_displacement_T(select_dofs(U, E, F)), epsv_times_h);
    }

protected:
    virtual double compute_distance(const VectorMax12d& x) const override;
    virtual VectorMax12d
    compute_distance_gradient(const VectorMax12d& x) const override;

    MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& x) const override;

    MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& x) const override;

    VectorMax2d compute_closest_point(const VectorMax12d& x) const override;

    MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& x) const override;

    VectorMax3d relative_displacement(const VectorMax12d& u) const override
    {
        return relative_displacement_T(u);
    }

    using FrictionConstraint::relative_displacement_matrix;

    MatrixMax<double, 3, 12> relative_displacement_matrix(
        const VectorMax2d& closest_point) const override;

    MatrixMax<double, 6, 12> relative_displacement_matrix_jacobian(
        const VectorMax2d& closest_point) const override;

private:
    template <typename T>
    VectorMax3<T> relative_displacement_T(const VectorMax12<T>& u) const
    {
        assert(u.size() == ndof());
        return edge_edge_relative_displacement(
            u.head(dim()), u.segment(dim(), dim()), u.segment(2 * dim(), dim()),
            u.tail(dim()), closest_point.cast<T>());
    }
};

///////////////////////////////////////////////////////////////////////////////

struct FaceVertexFrictionConstraint : FaceVertexCandidate, FrictionConstraint {
    FaceVertexFrictionConstraint(long face_index, long vertex_index);
    FaceVertexFrictionConstraint(const FaceVertexCandidate& constraint);
    FaceVertexFrictionConstraint(const FaceVertexConstraint& constraint);
    FaceVertexFrictionConstraint(
        const FaceVertexConstraint& constraint,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const double barrier_stiffness)
        : FaceVertexFrictionConstraint(constraint)
    {
        FrictionConstraint::init(
            V, E, F, dhat, barrier_stiffness, constraint.minimum_distance);
    }

    int num_vertices() const override { return 4; }
    std::array<long, 4> vertex_indices(
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
            relative_displacement_T(select_dofs(U, E, F)), epsv_times_h);
    }

protected:
    virtual double compute_distance(const VectorMax12d& x) const override;
    virtual VectorMax12d
    compute_distance_gradient(const VectorMax12d& x) const override;

    MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& x) const override;

    MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& x) const override;

    VectorMax2d compute_closest_point(const VectorMax12d& x) const override;

    MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& x) const override;

    VectorMax3d relative_displacement(const VectorMax12d& u) const override
    {
        return relative_displacement_T(u);
    }

    using FrictionConstraint::relative_displacement_matrix;

    MatrixMax<double, 3, 12> relative_displacement_matrix(
        const VectorMax2d& closest_point) const override;

    MatrixMax<double, 6, 12> relative_displacement_matrix_jacobian(
        const VectorMax2d& closest_point) const override;

private:
    template <typename T>
    VectorMax3<T> relative_displacement_T(const VectorMax12<T>& u) const
    {
        assert(u.size() == ndof());
        return point_triangle_relative_displacement(
            u.head(dim()), u.segment(dim(), dim()), u.segment(2 * dim(), dim()),
            u.tail(dim()), closest_point.cast<T>());
    }
};

///////////////////////////////////////////////////////////////////////////////

struct FrictionConstraints {
    std::vector<VertexVertexFrictionConstraint> vv_constraints;
    std::vector<EdgeVertexFrictionConstraint> ev_constraints;
    std::vector<EdgeEdgeFrictionConstraint> ee_constraints;
    std::vector<FaceVertexFrictionConstraint> fv_constraints;

    size_t size() const;

    bool empty() const;

    void clear();

    FrictionConstraint& operator[](size_t idx);
    const FrictionConstraint& operator[](size_t idx) const;
};

} // namespace ipc
