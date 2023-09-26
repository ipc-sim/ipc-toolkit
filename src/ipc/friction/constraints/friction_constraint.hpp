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

    /// @brief Compute the friction dissapative potential.
    /// @param velocities Velocities of the vertices (rowwise)
    /// @param edges Edges of the mesh
    /// @param faces Faces of the mesh
    /// @param epsv Smooth friction mollifier parameter \f$\epsilon_v\f$.
    /// @return The friction dissapative potential.
    double compute_potential(
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double epsv) const;

    /// @brief Compute the friction dissapative potential gradient wrt velocities.
    /// @param velocities Velocities of the vertices (rowwise)
    /// @param edges Edges of the mesh
    /// @param faces Faces of the mesh
    /// @param epsv Smooth friction mollifier parameter \f$\epsilon_v\f$.
    /// @return Gradient of the friction dissapative potential wrt velocities
    VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double epsv) const;

    /// @brief Compute the friction dissapative potential hessian wrt velocities.
    /// @param velocities Velocities of the vertices (rowwise)
    /// @param edges Edges of the mesh
    /// @param faces Faces of the mesh
    /// @param epsv Smooth friction mollifier parameter \f$\epsilon_v\f$.
    /// @param project_hessian_to_psd Project the hessian to PSD
    /// @return Hessian of the friction dissapative potential wrt velocities
    MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double epsv,
        const bool project_hessian_to_psd) const;

    /// @brief Compute the friction force.
    /// @param rest_positions Rest positions of the vertices (rowwise)
    /// @param lagged_displacements Previous displacements of the vertices (rowwise)
    /// @param velocities Current displacements of the vertices (rowwise)
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param epsv Smooth friction mollifier parameter \f$\epsilon_v\f$.
    /// @param dmin Minimum distance
    /// @param no_mu Whether to not multiply by mu
    /// @return Friction force
    VectorMax12d compute_force(
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& lagged_displacements,
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv,
        const double dmin = 0,
        const bool no_mu = false) const; //< whether to not multiply by mu

    /// @brief Variable to differentiate the friction force with respect to.
    enum class DiffWRT {
        REST_POSITIONS,       ///< Differentiate w.r.t. rest positions
        LAGGED_DISPLACEMENTS, ///< Differentiate w.r.t. lagged displacements
        VELOCITIES            ///< Differentiate w.r.t. current velocities
    };

    /// @brief Compute the friction force Jacobian.
    /// @param rest_positions Rest positions of the vertices (rowwise)
    /// @param lagged_displacements Previous displacements of the vertices (rowwise)
    /// @param velocities Current displacements of the vertices (rowwise)
    /// @param edges Collision mesh edges
    /// @param faces Collision mesh faces
    /// @param dhat Barrier activation distance
    /// @param barrier_stiffness Barrier stiffness
    /// @param epsv Smooth friction mollifier parameter \f$\epsilon_v\f$.
    /// @param wrt Variable to differentiate the friction force with respect to.
    /// @param dmin Minimum distance
    /// @return Friction force Jacobian
    MatrixMax12d compute_force_jacobian(
        const Eigen::MatrixXd& rest_positions,
        const Eigen::MatrixXd& lagged_displacements,
        const Eigen::MatrixXd& velocities,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const double dhat,
        const double barrier_stiffness,
        const double epsv,
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
    /// @return Relative velocity of the constraint.
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
