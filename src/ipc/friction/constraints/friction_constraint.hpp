#pragma once

#include <ipc/collisions/collision_constraints.hpp>
#include <ipc/friction/relative_velocity.hpp>
#include <ipc/friction/smooth_friction_mollifier.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <ipc/config.hpp>

namespace ipc {

class FrictionConstraint : virtual public CollisionStencil {
protected:
    /// @brief Initialize the constraint.
    /// @param positions Vertex positions(rowwise)
    /// @param edges Edges of the mesh
    /// @param faces Faces of the mesh
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param dmin Minimum distance
    void init(
        const Eigen::MatrixXd& positions,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double dmin);

public:
    virtual ~FrictionConstraint() { }

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

    /// @brief Compute the friction force.
    /// @param X Rest positions of the vertices (rowwise)
    /// @param U Current displacements of the vertices (rowwise)
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param epsv_times_h $\epsilon_vh$
    /// @param dmin Minimum distance
    /// @param no_mu Whether to not multiply by mu
    /// @return Friction force
    virtual VectorMax12d compute_force(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) const //< whether to not multiply by mu
    {
        return compute_force(
            X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U, edges, faces, dhat,
            barrier_stiffness, epsv_times_h, dmin, no_mu);
    }

    /// @brief Compute the friction force.
    /// @param X Rest positions of the vertices (rowwise)
    /// @param Ut Previous displacements of the vertices (rowwise)
    /// @param U Current displacements of the vertices (rowwise)
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param epsv_times_h $\epsilon_vh$
    /// @param dmin Minimum distance
    /// @param no_mu Whether to not multiply by mu
    /// @return Friction force
    virtual VectorMax12d compute_force(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const double dmin = 0,
        const bool no_mu = false) const; //< whether to not multiply by mu

    /// @brief Variable to differentiate the friction force with respect to.
    enum class DiffWRT {
        X,  ///< Differentiate wrt rest positions.
        Ut, ///< Differentiate wrt previous displacements.
        U   ///< Differentiate wrt current displacements.
    };

    /// @brief Compute the friction force Jacobian.
    /// @param X Rest positions of the vertices (rowwise)
    /// @param U Current displacements of the vertices (rowwise)
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param epsv_times_h $\epsilon_vh$
    /// @param wrt Variable to differentiate the friction force with respect to.
    /// @param dmin Minimum distance
    /// @return Friction force Jacobian
    virtual MatrixMax12d compute_force_jacobian(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& U,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv_times_h,
        const DiffWRT wrt,
        const double dmin = 0) const
    {
        return compute_force_jacobian(
            X, Eigen::MatrixXd::Zero(U.rows(), U.cols()), U, edges, faces, dhat,
            barrier_stiffness, epsv_times_h, wrt, dmin);
    }

    /// @brief Compute the friction force Jacobian.
    /// @param X Rest positions of the vertices (rowwise)
    /// @param Ut Previous displacements of the vertices (rowwise)
    /// @param U Current displacements of the vertices (rowwise)
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param epsv_times_h $\epsilon_vh$
    /// @param wrt Variable to differentiate the friction force with respect to.
    /// @param dmin Minimum distance
    /// @return Friction force Jacobian
    virtual MatrixMax12d compute_force_jacobian(
        const Eigen::MatrixXd& X,
        const Eigen::MatrixXd& Ut,
        const Eigen::MatrixXd& U,
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

    // -------------------------------------------------------------------------
    // Abstract methods
    // -------------------------------------------------------------------------

    /// @brief Compute the normal force magnitude.
    /// @param positions Constraint's vertex positions.
    /// @param dhat Barrier activiation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param dmin Minimum distance.
    /// @return Normal force magnitude.
    double compute_normal_force_magnitude(
        const VectorMax12d& positions,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const;

    /// @brief Compute the gradient of the normal force magnitude.
    /// @param positions Constraint's vertex positions.
    /// @param dhat Barrier activiation distance.
    /// @param barrier_stiffness Barrier stiffness.
    /// @param dmin Minimum distance.
    /// @return Gradient of the normal force magnitude wrt positions.
    VectorMax12d compute_normal_force_magnitude_gradient(
        const VectorMax12d& positions,
        const double dhat,
        const double barrier_stiffness,
        const double dmin = 0) const;

    /// @brief Compute the tangent basis of the constraint.
    /// @param positions Constraint's vertex positions.
    /// @return Tangent basis of the constraint.
    virtual MatrixMax<double, 3, 2>
    compute_tangent_basis(const VectorMax12d& positions) const = 0;

    /// @brief Compute the Jacobian of the tangent basis of the constraint.
    /// @param positions Constraint's vertex positions.
    /// @return Jacobian of the tangent basis of the constraint.
    virtual MatrixMax<double, 36, 2>
    compute_tangent_basis_jacobian(const VectorMax12d& positions) const = 0;

    /// @brief Compute the barycentric coordinates of the closest point.
    /// @param positions Constraint's vertex positions.
    /// @return Barycentric coordinates of the closest point.
    virtual VectorMax2d
    compute_closest_point(const VectorMax12d& positions) const = 0;

    /// @brief Compute the Jacobian of the barycentric coordinates of the closest point.
    /// @param positions Constraint's vertex positions.
    /// @return Jacobian of the barycentric coordinates of the closest point.
    virtual MatrixMax<double, 2, 12>
    compute_closest_point_jacobian(const VectorMax12d& positions) const = 0;

    /// @brief Compute the relative velocity of the constraint.
    /// @param velocities Constraint's vertex velocities.
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
