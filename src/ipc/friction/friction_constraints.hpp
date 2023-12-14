#pragma once

#include <ipc/friction/constraints/friction_constraint.hpp>
#include <ipc/friction/constraints/vertex_vertex.hpp>
#include <ipc/friction/constraints/edge_vertex.hpp>
#include <ipc/friction/constraints/edge_edge.hpp>
#include <ipc/friction/constraints/face_vertex.hpp>

#include <ipc/collision_mesh.hpp>
#include <ipc/collisions/collision_constraints.hpp>
#include <ipc/utils/eigen_ext.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

namespace ipc {

class FrictionConstraints {
public:
    using value_type = FrictionConstraint;

public:
    FrictionConstraints() { }

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const CollisionConstraints& contact_constraints,
        double dhat,
        double barrier_stiffness,
        double mu)
    {
        this->build(
            mesh, vertices, contact_constraints, dhat, barrier_stiffness,
            Eigen::VectorXd::Constant(vertices.rows(), mu));
    }

    void build(
        const CollisionMesh& mesh,
        const Eigen::MatrixXd& vertices,
        const CollisionConstraints& contact_constraints,
        const double dhat,
        const double barrier_stiffness,
        const Eigen::VectorXd& mus,
        const std::function<double(double, double)>& blend_mu =
            default_blend_mu);

    // ------------------------------------------------------------------------

    /// @brief Get the number of friction constraints.
    size_t size() const;

    /// @brief Get if the friction constraints are empty.
    bool empty() const;

    /// @brief Clear the friction constraints.
    void clear();

    /// @brief Get a reference to constriant idx.
    /// @param idx The index of the constraint.
    /// @return A reference to the constraint.
    FrictionConstraint& operator[](const size_t idx);

    /// @brief Get a const reference to constriant idx.
    /// @param idx The index of the constraint.
    /// @return A const reference to the constraint.
    const FrictionConstraint& operator[](const size_t idx) const;

    static double default_blend_mu(double mu0, double mu1)
    {
        // return mu0 * mu1;
        // return std::min(mu0, mu1);
        // return std::max(mu0, mu1);
        return (mu0 + mu1) / 2;
    }

public:
    std::vector<VertexVertexFrictionConstraint> vv_constraints;
    std::vector<EdgeVertexFrictionConstraint> ev_constraints;
    std::vector<EdgeEdgeFrictionConstraint> ee_constraints;
    std::vector<FaceVertexFrictionConstraint> fv_constraints;
};

} // namespace ipc
