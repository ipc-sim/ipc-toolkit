#pragma once

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct CollisionConstraint {
    virtual ~CollisionConstraint() {}

    virtual std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const = 0;

    virtual double compute_potential(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const = 0;

    virtual VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const = 0;

    virtual MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const bool project_hessian_to_psd) const = 0;

    double minimum_distance = 0;
};

///////////////////////////////////////////////////////////////////////////////

struct VertexVertexConstraint : VertexVertexCandidate, CollisionConstraint {
    VertexVertexConstraint(long vertex0_index, long vertex1_index);
    VertexVertexConstraint(const VertexVertexCandidate& candidate);

    std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex0_index, vertex1_index } };
    }

    double compute_potential(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const bool project_hessian_to_psd) const override;

    long multiplicity = 1;
};

///////////////////////////////////////////////////////////////////////////////

struct EdgeVertexConstraint : EdgeVertexCandidate, CollisionConstraint {
    EdgeVertexConstraint(long edge_index, long vertex_index);
    EdgeVertexConstraint(const EdgeVertexCandidate& candidate);

    std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, E(edge_index, 0), E(edge_index, 1) } };
    }

    double compute_potential(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const bool project_hessian_to_psd) const override;

    long multiplicity = 1;
};

///////////////////////////////////////////////////////////////////////////////

struct EdgeEdgeConstraint : EdgeEdgeCandidate, CollisionConstraint {
    EdgeEdgeConstraint(long edge0_index, long edge1_index, double eps_x);
    EdgeEdgeConstraint(const EdgeEdgeCandidate& candidate, double eps_x);

    std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { E(edge0_index, 0), E(edge0_index, 1), //
                   E(edge1_index, 0), E(edge1_index, 1) } };
    }

    double compute_potential(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const bool project_hessian_to_psd) const override;

    double eps_x;
};

///////////////////////////////////////////////////////////////////////////////

struct FaceVertexConstraint : FaceVertexCandidate, CollisionConstraint {
    FaceVertexConstraint(long face_index, long vertex_index);
    FaceVertexConstraint(const FaceVertexCandidate& candidate);

    std::vector<long> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, //
                   F(face_index, 0), F(face_index, 1), F(face_index, 2) } };
    }

    double compute_potential(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const override;

    MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const bool project_hessian_to_psd) const override;
};

///////////////////////////////////////////////////////////////////////////////

struct Constraints {
    std::vector<VertexVertexConstraint> vv_constraints;
    std::vector<EdgeVertexConstraint> ev_constraints;
    std::vector<EdgeEdgeConstraint> ee_constraints;
    std::vector<FaceVertexConstraint> fv_constraints;

    size_t size() const;

    size_t num_constraints() const;

    void clear();

    CollisionConstraint& operator[](size_t idx);
    const CollisionConstraint& operator[](size_t idx) const;
};

} // namespace ipc
