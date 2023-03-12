#pragma once

#include <ipc/collisions/collision_constraints.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <ipc/config.hpp>

namespace ipc {

class FrictionConstraint {
protected:
    /// @brief Initialize the constraint.
    /// @param V Positions of the vertices (rowwise)
    /// @param edges Edges of the mesh
    /// @param faces Faces of the mesh
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param dmin Minimum distance
    void init(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double dmin);

public:
    virtual ~FrictionConstraint() { }

    /// @brief Get the number of vertices in the friction constraint.
    virtual int num_vertices() const = 0;

    /// @brief Get the indices of the vertices in the friction constraint.
    virtual std::array<long, 4> vertex_ids(
        const Eigen::MatrixXi& edges, const Eigen::MatrixXi& faces) const = 0;

    /// @brief Compute the friction dissapative potential gradient wrt velocities.
    /// @param velocities Velocities of the vertices (rowwise)
    /// @param edges Edges of the mesh
    /// @param faces Faces of the mesh
    /// @param epsv_times_h $\epsilon_vh$
    /// @return Gradient of the friction dissapative potential wrt velocities
    virtual VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double epsv_times_h) const;

    /// @brief Compute the friction dissapative potential hessian wrt velocities.
    /// @param velocities Velocities of the vertices (rowwise)
    /// @param edges Edges of the mesh
    /// @param faces Faces of the mesh
    /// @param epsv_times_h $\epsilon_vh$
    /// @param project_hessian_to_psd Project the hessian to PSD
    /// @return Hessian of the friction dissapative potential wrt velocities
    virtual MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double epsv_times_h,
        const bool project_hessian_to_psd) const;

    /// @brief
    /// @param X
    /// @param velocities
    /// @param edges
    /// @param faces
    /// @param dhat
    /// @param barrier_stiffness
    /// @param epsv_times_h
    /// @param dmin
    /// @param no_mu
    /// @return
    virtual VectorMax12d compute_force(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) const //< whether to not multiply by mu
    {
        return compute_force(
            X, Eigen::MatrixXd::Zero(velocities.rows(), velocities.cols()),
            velocities, edges, faces, dhat, barrier_stiffness, epsv_times_h,
            dmin, no_mu);
    }

    /// @brief
    /// @param X
    /// @param Ut
    /// @param velocities
    /// @param edges
    /// @param faces
    /// @param dhat
    /// @param barrier_stiffness
    /// @param epsv_times_h
    /// @param dmin
    /// @param no_mu
    /// @return
    virtual VectorMax12d compute_force(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) const; //< whether to not multiply by mu

    /// @brief Variable to differentiate the friction force with respect to.
    enum class DiffWRT {
        X,  ///< Differentiate wrt rest V.
        Ut, ///< Differentiate wrt previous V.
        U   ///< Differentiate wrt current V.
    };

    /// @brief
    /// @param X
    /// @param velocities
    /// @param edges
    /// @param faces
    /// @param dhat
    /// @param barrier_stiffness
    /// @param epsv_times_h
    /// @param wrt
    /// @param dmin
    /// @return
    virtual MatrixMax12d compute_force_jacobian(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const DiffWRT wrt,
        const double dmin = 0) const
    {
        return compute_force_jacobian(
            X, Eigen::MatrixXd::Zero(velocities.rows(), velocities.cols()),
            velocities, edges, faces, dhat, barrier_stiffness, epsv_times_h,
            wrt, dmin);
    }

    /// @brief
    /// @param X
    /// @param Ut
    /// @param velocities
    /// @param edges
    /// @param faces
    /// @param dhat
    /// @param barrier_stiffness
    /// @param epsv_times_h
    /// @param wrt
    /// @param dmin
    /// @return
    virtual MatrixMax12d compute_force_jacobian(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const DiffWRT wrt,
        const double dmin = 0) const;

protected:
    /// @brief Get the dimension of the constraint.
    int dim() const { return tangent_basis.rows(); }

    /// @brief Get the number of degrees of freedom for the constraint.
    int ndof() const { return dim() * num_vertices(); };

    /// @brief Select this constraint's DOF from the full matrix of DOF.
    /// @tparam T Type of the DOF
    /// @param X Full matrix of DOF (rowwise).
    /// @param edges Edges of the mesh.
    /// @param faces Faces of the mesh.
    /// @return This constraint's DOF.
    template <typename T>
    VectorMax12<T> select_dof(
        const MatrixX<T>& DOF,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces) const
    {
        VectorMax12<T> dof(ndof());
        const std::array<long, 4> idx = vertex_ids(edges, faces);
        for (int i = 0; i < num_vertices(); i++) {
            dof.segment(i * dim(), dim()) = DOF.row(idx[i]);
        }
        return dof;
    }

    // -------------------------------------------------------------------------
    // Abstract methods
    // -------------------------------------------------------------------------

    /// @brief Compute the distance of the constraint.
    /// @param V Vertex V for the constraint's vertices (stacked into a vector).
    /// @return Distance of the constraint.
    virtual double compute_distance(const VectorMax12d& V) const = 0;

    /// @brief Compute the gradient of the distance of the constraint.
    /// @param V Vertex V for the constraint's vertices (stacked into a vector).
    /// @return Gradient of the distance of the constraint.
    virtual VectorMax12d
    compute_distance_gradient(const VectorMax12d& V) const = 0;

    /// @brief Compute the normal force magnitude.
    /// @param V Vertex V for the constraint's vertices (stacked into a vector).
    /// @param dhat Barrier activiation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param dmin Minimum distance.
    /// @return Normal force magnitude.
    virtual double compute_normal_force_magnitude(
        const VectorMax12d& V,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const;

    /// @brief Compute the gradient of the normal force magnitude.
    /// @param V Vertex V for the constraint's vertices (stacked into a vector).
    /// @param dhat Barrier activiation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param dmin Minimum distance.
    /// @return Gradient of the normal force magnitude wrt V.
    virtual VectorMax12d compute_normal_force_magnitude_gradient(
        const VectorMax12d& V,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const;

    /// @brief Compute the tangent basis of the constraint.
    /// @param V Vertex V for the constraint's vertices (stacked into a vector).
    /// @return Tangent basis of the constraint.
    virtual MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& V) const = 0;

    /// @brief Compute the Jacobian of the tangent basis of the constraint.
    /// @param V Vertex V for the constraint's vertices (stacked into a vector).
    /// @return Jacobian of the tangent basis of the constraint.
    virtual MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& V) const = 0;

    /// @brief Compute the barycentric coordinates of the closest point.
    /// @param V Vertex V for the constraint's vertices (stacked into a vector).
    /// @return Barycentric coordinates of the closest point.
    virtual VectorMax2d compute_closest_point(const VectorMax12d& V) const = 0;

    /// @brief Compute the Jacobian of the barycentric coordinates of the closest point.
    /// @param V Vertex V for the constraint's vertices (stacked into a vector).
    /// @return Jacobian of the barycentric coordinates of the closest point.
    virtual MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& V) const = 0;

    /// @brief Compute the relative velocity of the constraint.
    /// @param velocities Vertex velocities for the constraint's vertices (stacked into a vector).
    virtual VectorMax3d
    relative_velocity(const VectorMax12d& velocities) const = 0;

    /// @brief Construct the premultiplier matrix for the relative velocity.
    /// @note Uses the cached closest point.
    /// @return A matrix M such that `relative_velocity = M * velocities`.
    virtual MatrixMax<double, 3, 12> relative_velocity_matrix() const
    {
        return relative_velocity_matrix(closest_point);
    }

    /// @brief Construct the premultiplier matrix for the relative velocity.
    /// @param closest_point Barycentric coordinates of the closest point.
    /// @return A matrix M such that `relative_velocity = M * velocities`.
    virtual MatrixMax<double, 3, 12>
    relative_velocity_matrix(const VectorMax2d& closest_point) const = 0;

    /// @brief Construct the Jacobian of the relative velocity premultiplier wrt the closest points.
    /// @param closest_point Barycentric coordinates of the closest point.
    /// @return Jacobian of the relative velocity premultiplier wrt the closest points.
    virtual MatrixMax<double, 6, 12> relative_velocity_matrix_jacobian(
        const VectorMax2d& closest_point) const = 0;

    // -------------------------------------------------------------------------
    // Utility methods
    // -------------------------------------------------------------------------

    /// @brief Compute the friction potential from the relative velocity.
    /// @param relative_velocity Relative velocity of the constraint.
    /// @param epsv_times_h Friction mollifier parameter.
    /// @return Friction potential.
    template <typename DerivedRelUi, typename T = typename DerivedRelUi::Scalar>
    T compute_potential_common(
        const Eigen::MatrixBase<DerivedRelUi>& relative_velocity,
        const double epsv_times_h) const
    {
        // The relative velocity in the tangential space
        const VectorMax2d tangent_rel_vel =
            tangent_basis.transpose().cast<T>() * relative_velocity;
        return weight * mu * normal_force_magnitude
            * f0_SF(tangent_rel_vel.norm(), epsv_times_h);
    }

public:
    /// @brief Contact force magnitude
    double normal_force_magnitude;

    /// @brief Coefficient of friction
    double mu;

    /// @brief Weight
    double weight = 1;

    /// @brief Gradient of weight with respect to all DOF
    Eigen::SparseVector<double> weight_gradient;

    /// @brief Barycentric coordinates of the closest point(s)
    VectorMax2d closest_point;

    /// @brief Tangent basis of the contact (max size 3×2)
    MatrixMax<double, 3, 2> tangent_basis;
};

} // namespace ipc
