#pragma once

#include <array>

#include <Eigen/Core>

#include <ipc/broad_phase/collision_candidate.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct CollisionConstraint {
public:
    virtual ~CollisionConstraint() { }

    virtual int num_vertices() const = 0;
    virtual std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const = 0;

    virtual double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const = 0;

    virtual double compute_potential(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const;

    virtual VectorMax12d compute_potential_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat) const;

    virtual MatrixMax12d compute_potential_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F,
        const double dhat,
        const bool project_hessian_to_psd) const;

    double minimum_distance = 0;
};

///////////////////////////////////////////////////////////////////////////////

struct VertexVertexConstraint : VertexVertexCandidate, CollisionConstraint {
    VertexVertexConstraint(long vertex0_index, long vertex1_index);
    VertexVertexConstraint(const VertexVertexCandidate& candidate);

    int num_vertices() const override { return 2; };
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex0_index, vertex1_index, -1, -1 } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

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

    template <typename H>
    friend H AbslHashValue(H h, const VertexVertexConstraint& vv)
    {
        long min_vi = std::min(vv.vertex0_index, vv.vertex1_index);
        long max_vi = std::max(vv.vertex0_index, vv.vertex1_index);
        return H::combine(std::move(h), min_vi, max_vi);
    }

    int multiplicity = 1;
};

///////////////////////////////////////////////////////////////////////////////

struct EdgeVertexConstraint : EdgeVertexCandidate, CollisionConstraint {
    EdgeVertexConstraint(long edge_index, long vertex_index);
    EdgeVertexConstraint(const EdgeVertexCandidate& candidate);

    int num_vertices() const override { return 3; };
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, E(edge_index, 0), E(edge_index, 1), -1 } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

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

    template <typename H>
    friend H AbslHashValue(H h, const EdgeVertexConstraint& ev)
    {
        return H::combine(std::move(h), ev.edge_index, ev.vertex_index);
    }

    int multiplicity = 1;
};

///////////////////////////////////////////////////////////////////////////////

struct EdgeEdgeConstraint : EdgeEdgeCandidate, CollisionConstraint {
    EdgeEdgeConstraint(long edge0_index, long edge1_index, double eps_x);
    EdgeEdgeConstraint(const EdgeEdgeCandidate& candidate, double eps_x);

    int num_vertices() const override { return 4; };
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { E(edge0_index, 0), E(edge0_index, 1), //
                   E(edge1_index, 0), E(edge1_index, 1) } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

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

    template <typename H>
    friend H AbslHashValue(H h, const EdgeEdgeConstraint& ee)
    {
        long min_ei = std::min(ee.edge0_index, ee.edge1_index);
        long max_ei = std::max(ee.edge0_index, ee.edge1_index);
        return H::combine(std::move(h), min_ei, max_ei);
    }

    double eps_x;
};

///////////////////////////////////////////////////////////////////////////////

struct FaceVertexConstraint : FaceVertexCandidate, CollisionConstraint {
    FaceVertexConstraint(long face_index, long vertex_index);
    FaceVertexConstraint(const FaceVertexCandidate& candidate);

    int num_vertices() const override { return 4; };
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, //
                   F(face_index, 0), F(face_index, 1), F(face_index, 2) } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    template <typename H>
    friend H AbslHashValue(H h, const FaceVertexConstraint& fv)
    {
        return H::combine(std::move(h), fv.face_index, fv.vertex_index);
    }
};

///////////////////////////////////////////////////////////////////////////////

struct PlaneVertexConstraint : CollisionConstraint {
    PlaneVertexConstraint(
        const VectorMax3d& plane_origin,
        const VectorMax3d& plane_normal,
        const long vertex_index);

    int num_vertices() const override { return 1; };
    std::array<long, 4> vertex_indices(
        const Eigen::MatrixXi& E, const Eigen::MatrixXi& F) const override
    {
        return { { vertex_index, -1, -1, -1 } };
    }

    double compute_distance(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax12d compute_distance_gradient(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    MatrixMax12d compute_distance_hessian(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F) const override;

    VectorMax3d plane_origin;
    VectorMax3d plane_normal;
    long vertex_index;
};

///////////////////////////////////////////////////////////////////////////////

struct Constraints {
    std::vector<VertexVertexConstraint> vv_constraints;
    std::vector<EdgeVertexConstraint> ev_constraints;
    std::vector<EdgeEdgeConstraint> ee_constraints;
    std::vector<FaceVertexConstraint> fv_constraints;
    std::vector<PlaneVertexConstraint> pv_constraints;

    size_t size() const;

    size_t num_constraints() const;

    bool empty() const;

    void clear();

    CollisionConstraint& operator[](size_t idx);
    const CollisionConstraint& operator[](size_t idx) const;
};

} // namespace ipc

namespace std {
template <> struct hash<ipc::VertexVertexConstraint> {
    size_t operator()(ipc::VertexVertexConstraint const& vv) const noexcept;
};
template <> struct hash<ipc::EdgeVertexConstraint> {
    size_t operator()(ipc::EdgeVertexConstraint const& ev) const noexcept;
};
template <> struct hash<ipc::EdgeEdgeConstraint> {
    size_t operator()(ipc::EdgeEdgeConstraint const& ee) const noexcept;
};
template <> struct hash<ipc::FaceVertexConstraint> {
    size_t operator()(ipc::FaceVertexConstraint const& fv) const noexcept;
};
} // namespace std
