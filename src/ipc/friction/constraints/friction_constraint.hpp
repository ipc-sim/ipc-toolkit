#pragma once

#include <ipc/collisions/constraints.hpp>
#include <ipc/friction/relative_velocity.hpp>
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

    /// @brief Compute the friction dissapative potential gradient wrt U.
    /// @param U Velocities of the vertices (rowwise)
    /// @param E Edges of the mesh
    /// @param F Faces of the mesh
    /// @param epsv_times_h $\epsilon_vh$
    /// @return Gradient of the friction dissapative potential wrt U
    virtual VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double epsv_times_h) const;

    /// @brief Compute the friction dissapative potential hessian wrt U.
    /// @param U Velocities of the vertices (rowwise)
    /// @param E Edges of the mesh
    /// @param F Faces of the mesh
    /// @param epsv_times_h $\epsilon_vh$
    /// @param project_hessian_to_psd Project the hessian to PSD
    /// @return Hessian of the friction dissapative potential wrt U
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

    virtual VectorMax3d relative_velocity(const VectorMax12d& u) const = 0;

    virtual MatrixMax<double, 3, 12> relative_velocity_matrix() const
    {
        return relative_velocity_matrix(closest_point);
    }

    virtual MatrixMax<double, 3, 12>
    relative_velocity_matrix(const VectorMax2d& closest_point) const = 0;

    virtual MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        const VectorMax2d& closest_point) const = 0;

    template <typename DerivedRelUi, typename T = typename DerivedRelUi::Scalar>
    T compute_potential_common(
        const Eigen::MatrixBase<DerivedRelUi>& rel_ui,
        const double epsv_times_h) const
    {
        // u is the relative velocity in the tangential space
        const VectorMax2d u = tangent_basis.transpose().cast<T>() * rel_ui;
        return weight * mu * normal_force_magnitude
            * f0_SF(u.norm(), epsv_times_h);
    }
};

} // namespace ipc
