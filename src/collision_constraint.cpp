#include <ipc/collision_constraint.hpp>

#include <ipc/barrier/barrier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

VertexVertexConstraint::VertexVertexConstraint(
    long vertex0_index, long vertex1_index)
    : VertexVertexCandidate(vertex0_index, vertex1_index)
{
}

VertexVertexConstraint::VertexVertexConstraint(
    const VertexVertexCandidate& candidate)
    : VertexVertexCandidate(candidate)
{
}

double VertexVertexConstraint::compute_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin) const
{
    double dhat_squared = dhat * dhat;
    double distance_sqr =
        point_point_distance(V.row(vertex0_index), V.row(vertex1_index));
    return multiplicity
        * barrier(distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
}

Eigen::VectorX12d VertexVertexConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin) const
{
    double dhat_squared = dhat * dhat;

    // ∇[m * b(d(x))] = m * b'(d(x)) * ∇d(x)
    const auto& p0 = V.row(vertex0_index);
    const auto& p1 = V.row(vertex1_index);

    Eigen::VectorX6d local_grad;
    point_point_distance_gradient(p0, p1, local_grad);

    double distance_sqr = point_point_distance(p0, p1);
    local_grad *= barrier_gradient(
        distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);

    return multiplicity * local_grad;
}

Eigen::MatrixXX12d VertexVertexConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin,
    const bool project_to_psd) const
{
    double dhat_squared = dhat * dhat;

    // ∇²[m * b(d(x))] = m * ∇(b'(d(x)) * ∇d(x))
    //                 = m * [b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)]
    const auto& p0 = V.row(vertex0_index);
    const auto& p1 = V.row(vertex1_index);

    double distance_sqr = point_point_distance(p0, p1);
    Eigen::VectorX6d local_grad;
    point_point_distance_gradient(p0, p1, local_grad);
    Eigen::MatrixXX6d local_hess;
    point_point_distance_hessian(p0, p1, local_hess);

    local_hess *= barrier_gradient(
        distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
    local_hess +=
        barrier_hessian(
            distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared)
        * local_grad * local_grad.transpose();

    local_hess *= multiplicity;

    if (project_to_psd) {
        local_hess = Eigen::project_to_psd(local_hess);
    }

    return local_hess;
}

///////////////////////////////////////////////////////////////////////////////

EdgeVertexConstraint::EdgeVertexConstraint(long edge_index, long vertex_index)
    : EdgeVertexCandidate(edge_index, vertex_index)
{
}

EdgeVertexConstraint::EdgeVertexConstraint(const EdgeVertexCandidate& candidate)
    : EdgeVertexCandidate(candidate)
{
}

double EdgeVertexConstraint::compute_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin) const
{
    double dhat_squared = dhat * dhat;
    // The distance type is known because of construct_constraint_set()
    double distance_sqr = point_edge_distance(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)),
        PointEdgeDistanceType::P_E);
    return multiplicity
        * barrier(distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
}

Eigen::VectorX12d EdgeVertexConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin) const
{
    double dhat_squared = dhat * dhat;

    // ∇[m * b(d(x))] = m * b'(d(x)) * ∇d(x)
    const auto& p = V.row(vertex_index);
    const auto& e0 = V.row(E(edge_index, 0));
    const auto& e1 = V.row(E(edge_index, 1));

    Eigen::VectorX9d local_grad;
    point_edge_distance_gradient(
        p, e0, e1, PointEdgeDistanceType::P_E, local_grad);

    double distance_sqr =
        point_edge_distance(p, e0, e1, PointEdgeDistanceType::P_E);
    local_grad *= barrier_gradient(
        distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);

    return multiplicity * local_grad;
}

Eigen::MatrixXX12d EdgeVertexConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin,
    const bool project_to_psd) const
{
    double dhat_squared = dhat * dhat;

    // ∇²[m * b(d(x))] = m * ∇(b'(d(x)) * ∇d(x))
    //                 = m * [b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)]
    const auto& p = V.row(vertex_index);
    const auto& e0 = V.row(E(edge_index, 0));
    const auto& e1 = V.row(E(edge_index, 1));

    double distance_sqr =
        point_edge_distance(p, e0, e1, PointEdgeDistanceType::P_E);
    Eigen::VectorX9d local_grad;
    point_edge_distance_gradient(
        p, e0, e1, PointEdgeDistanceType::P_E, local_grad);
    Eigen::MatrixXX12d local_hess;
    point_edge_distance_hessian(
        p, e0, e1, PointEdgeDistanceType::P_E, local_hess);

    local_hess *= barrier_gradient(
        distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
    local_hess +=
        barrier_hessian(
            distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared)
        * local_grad * local_grad.transpose();

    local_hess *= multiplicity;

    if (project_to_psd) {
        local_hess = Eigen::project_to_psd(local_hess);
    }

    return local_hess;
}

///////////////////////////////////////////////////////////////////////////////

EdgeEdgeConstraint::EdgeEdgeConstraint(
    long edge0_index, long edge1_index, double eps_x)
    : EdgeEdgeCandidate(edge0_index, edge1_index)
    , eps_x(eps_x)
{
}

EdgeEdgeConstraint::EdgeEdgeConstraint(
    const EdgeEdgeCandidate& candidate, double eps_x)
    : EdgeEdgeCandidate(candidate)
    , eps_x(eps_x)
{
}

double EdgeEdgeConstraint::compute_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin) const
{
    double dhat_squared = dhat * dhat;

    const auto& ea0 = V.row(E(edge0_index, 0));
    const auto& ea1 = V.row(E(edge0_index, 1));
    const auto& eb0 = V.row(E(edge1_index, 0));
    const auto& eb1 = V.row(E(edge1_index, 1));

    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1);
    return edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x)
        * barrier(distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
}

Eigen::VectorX12d EdgeEdgeConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin) const
{
    double dhat_squared = dhat * dhat;

    // ∇[m(x) * b(d(x))] = (∇m(x)) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)
    const auto& ea0 = V.row(E(edge0_index, 0));
    const auto& ea1 = V.row(E(edge0_index, 1));
    const auto& eb0 = V.row(E(edge1_index, 0));
    const auto& eb1 = V.row(E(edge1_index, 1));

    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    EdgeEdgeDistanceType dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
    Eigen::VectorX12d local_distance_grad;
    edge_edge_distance_gradient(ea0, ea1, eb0, eb1, dtype, local_distance_grad);

    double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    Eigen::VectorX12d local_mollifier_grad;
    edge_edge_mollifier_gradient(
        ea0, ea1, eb0, eb1, eps_x, local_mollifier_grad);

    return local_mollifier_grad
        * barrier(distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared)
        + mollifier
        * barrier_gradient(
              distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared)
        * local_distance_grad;
}

Eigen::MatrixXX12d EdgeEdgeConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin,
    const bool project_to_psd) const
{
    double dhat_squared = dhat * dhat;

    // ∇²[m(x) * b(d(x))] = ∇[∇m(x) * b(d(x)) + m(x) * b'(d(x)) * ∇d(x)]
    //                    = ∇²m(x) * b(d(x)) + b'(d(x)) * ∇d(x) * ∇m(x)ᵀ
    //                      + ∇m(x) * b'(d(x)) * ∇d(x))ᵀ
    //                      + m(x) * b"(d(x)) * ∇d(x) * ∇d(x)ᵀ
    //                      + m(x) * b'(d(x)) * ∇²d(x)
    const auto& ea0 = V.row(E(edge0_index, 0));
    const auto& ea1 = V.row(E(edge0_index, 1));
    const auto& eb0 = V.row(E(edge1_index, 0));
    const auto& eb1 = V.row(E(edge1_index, 1));

    // Compute distance derivatives
    // The distance type is unknown because of mollified PP and PE
    // constraints where also added as EE constraints.
    EdgeEdgeDistanceType dtype = edge_edge_distance_type(ea0, ea1, eb0, eb1);
    double distance_sqr = edge_edge_distance(ea0, ea1, eb0, eb1, dtype);
    Eigen::VectorX12d distance_grad;
    edge_edge_distance_gradient(ea0, ea1, eb0, eb1, dtype, distance_grad);
    Eigen::MatrixXX12d distance_hess;
    edge_edge_distance_hessian(ea0, ea1, eb0, eb1, dtype, distance_hess);

    // Compute mollifier derivatives
    double mollifier = edge_edge_mollifier(ea0, ea1, eb0, eb1, eps_x);
    Eigen::VectorX12d mollifier_grad;
    edge_edge_mollifier_gradient(ea0, ea1, eb0, eb1, eps_x, mollifier_grad);
    Eigen::MatrixXX12d mollifier_hess;
    edge_edge_mollifier_hessian(ea0, ea1, eb0, eb1, eps_x, mollifier_hess);

    // Compute_barrier_derivatives
    double b =
        barrier(distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
    double grad_b = barrier_gradient(
        distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
    double hess_b = barrier_hessian(
        distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);

    Eigen::MatrixXX12d local_hess = mollifier_hess * b
        + grad_b
            * (distance_grad * mollifier_grad.transpose()
               + mollifier_grad * distance_grad.transpose())
        + mollifier
            * (hess_b * distance_grad * distance_grad.transpose()
               + grad_b * distance_hess);

    if (project_to_psd) {
        local_hess = Eigen::project_to_psd(local_hess);
    }

    return local_hess;
}

///////////////////////////////////////////////////////////////////////////////

FaceVertexConstraint::FaceVertexConstraint(long face_index, long vertex_index)
    : FaceVertexCandidate(face_index, vertex_index)
{
}

FaceVertexConstraint::FaceVertexConstraint(const FaceVertexCandidate& candidate)
    : FaceVertexCandidate(candidate)
{
}

double FaceVertexConstraint::compute_potential(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin) const
{
    double dhat_squared = dhat * dhat;

    const auto& p = V.row(vertex_index);
    const auto& t0 = V.row(F(face_index, 0));
    const auto& t1 = V.row(F(face_index, 1));
    const auto& t2 = V.row(F(face_index, 2));

    // The distance type is known because of construct_constraint_set()
    double distance_sqr =
        point_triangle_distance(p, t0, t1, t2, PointTriangleDistanceType::P_T);
    return barrier(distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
}

Eigen::VectorX12d FaceVertexConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin) const
{
    double dhat_squared = dhat * dhat;

    // ∇b(d(x)) = b'(d(x)) * ∇d(x)
    const auto& p = V.row(vertex_index);
    const auto& t0 = V.row(F(face_index, 0));
    const auto& t1 = V.row(F(face_index, 1));
    const auto& t2 = V.row(F(face_index, 2));

    Eigen::VectorX12d local_grad;
    point_triangle_distance_gradient(
        p, t0, t1, t2, PointTriangleDistanceType::P_T, local_grad);

    double distance_sqr =
        point_triangle_distance(p, t0, t1, t2, PointTriangleDistanceType::P_T);

    return local_grad
        * barrier_gradient(
               distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
}

Eigen::MatrixXX12d FaceVertexConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double dmin,
    const bool project_to_psd) const
{
    double dhat_squared = dhat * dhat;

    // ∇²b(d(x)) = ∇(b'(d(x)) * ∇d(x))
    //           = b"(d(x)) * ∇d(x) * ∇d(x)ᵀ + b'(d(x)) * ∇²d(x)
    const auto& p = V.row(vertex_index);
    const auto& t0 = V.row(F(face_index, 0));
    const auto& t1 = V.row(F(face_index, 1));
    const auto& t2 = V.row(F(face_index, 2));

    double distance_sqr =
        point_triangle_distance(p, t0, t1, t2, PointTriangleDistanceType::P_T);
    Eigen::VectorX12d local_grad;
    point_triangle_distance_gradient(
        p, t0, t1, t2, PointTriangleDistanceType::P_T, local_grad);
    Eigen::MatrixXX12d local_hess;
    point_triangle_distance_hessian(
        p, t0, t1, t2, PointTriangleDistanceType::P_T, local_hess);

    local_hess *= barrier_gradient(
        distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared);
    local_hess +=
        barrier_hessian(
            distance_sqr - dmin * dmin, 2 * dmin * dhat + dhat_squared)
        * local_grad * local_grad.transpose();

    if (project_to_psd) {
        local_hess = Eigen::project_to_psd(local_hess);
    }

    return local_hess;
}

///////////////////////////////////////////////////////////////////////////////

size_t Constraints::size() const
{
    return vv_constraints.size() + ev_constraints.size() + ee_constraints.size()
        + fv_constraints.size();
}

size_t Constraints::num_constraints() const
{
    size_t num_constraints = 0;
    for (const auto& vv_constraint : vv_constraints) {
        num_constraints += vv_constraint.multiplicity;
    }
    for (const auto& ev_constraint : ev_constraints) {
        num_constraints += ev_constraint.multiplicity;
    }
    num_constraints += ee_constraints.size() + fv_constraints.size();
    return num_constraints;
}

void Constraints::clear()
{
    vv_constraints.clear();
    ev_constraints.clear();
    ee_constraints.clear();
    fv_constraints.clear();
}

CollisionConstraint& Constraints::operator[](size_t idx)
{
    if (idx < vv_constraints.size()) {
        return vv_constraints[idx];
    }
    idx -= vv_constraints.size();
    if (idx < ev_constraints.size()) {
        return ev_constraints[idx];
    }
    idx -= ev_constraints.size();
    if (idx < ee_constraints.size()) {
        return ee_constraints[idx];
    }
    idx -= ee_constraints.size();
    if (idx < fv_constraints.size()) {
        return fv_constraints[idx];
    }
    assert(false);
    throw "Invalid friction constraint index!";
}

const CollisionConstraint& Constraints::operator[](size_t idx) const
{
    if (idx < vv_constraints.size()) {
        return vv_constraints[idx];
    }
    idx -= vv_constraints.size();
    if (idx < ev_constraints.size()) {
        return ev_constraints[idx];
    }
    idx -= ev_constraints.size();
    if (idx < ee_constraints.size()) {
        return ee_constraints[idx];
    }
    idx -= ee_constraints.size();
    if (idx < fv_constraints.size()) {
        return fv_constraints[idx];
    }
    assert(false);
    throw "Invalid friction constraint index!";
}

} // namespace ipc
