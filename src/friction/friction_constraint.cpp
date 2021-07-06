#include <ipc/friction/friction_constraint.hpp>

#include <stdexcept> // std::out_of_range

namespace ipc {

VectorMax2d FrictionConstraint::compute_potential_gradient_common(
    const VectorMax3d& relative_displacement, double epsv_times_h) const
{
    VectorMax2d tangent_relative_displacement =
        tangent_basis.transpose() * relative_displacement;

    double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
        tangent_relative_displacement.squaredNorm(), epsv_times_h);

    tangent_relative_displacement *=
        f1_div_rel_disp_norm * mu * normal_force_magnitude;

    return tangent_relative_displacement;
}

MatrixMax12d FrictionConstraint::compute_potential_hessian_common(
    const VectorMax3d& relative_displacement,
    const MatrixMax<double, 2, 12>& TT,
    const double epsv_times_h,
    bool project_hessian_to_psd,
    const int multiplicity) const
{
    int dim = relative_displacement.size();
    assert(dim == 2 || dim == 3);

    double epsv_times_h_squared = epsv_times_h * epsv_times_h;

    VectorMax2d tangent_relative_displacement =
        tangent_basis.transpose() * relative_displacement;

    double tangent_relative_displacement_sqnorm =
        tangent_relative_displacement.squaredNorm();

    double f1_div_rel_disp_norm = f1_SF_div_relative_displacement_norm(
        tangent_relative_displacement_sqnorm, epsv_times_h);
    double f2_term = f2_SF(tangent_relative_displacement_sqnorm, epsv_times_h);

    MatrixMax12d local_hess;

    double scale = multiplicity * mu * normal_force_magnitude;
    if (tangent_relative_displacement_sqnorm >= epsv_times_h_squared) {
        // no SPD projection needed
        VectorMax2d ubar(dim - 1);
        if (dim == 2) {
            ubar[0] = tangent_relative_displacement[0];
        } else {
            ubar[0] = -tangent_relative_displacement[1];
            ubar[1] = tangent_relative_displacement[0];
        }
        local_hess = (TT.transpose()
                      * ((scale * f1_div_rel_disp_norm
                          / tangent_relative_displacement_sqnorm)
                         * ubar))
            * (ubar.transpose() * TT);
    } else {
        double tangent_relative_displacement_norm =
            sqrt(tangent_relative_displacement_sqnorm);
        if (tangent_relative_displacement_norm == 0) {
            // no SPD projection needed
            local_hess = ((scale * f1_div_rel_disp_norm) * TT.transpose()) * TT;
        } else {
            // only need to project the inner 2x2 matrix to SPD
            MatrixMax2d inner_hess =
                ((f2_term / tangent_relative_displacement_norm)
                 * tangent_relative_displacement)
                * tangent_relative_displacement.transpose();
            inner_hess.diagonal().array() += f1_div_rel_disp_norm;
            if (project_hessian_to_psd) {
                inner_hess = project_to_psd(inner_hess);
            }
            inner_hess *= scale;

            // tensor product:
            local_hess = TT.transpose() * inner_hess * TT;
        }
    }

    return local_hess;
}

///////////////////////////////////////////////////////////////////////////////

VertexVertexFrictionConstraint::VertexVertexFrictionConstraint(
    long vertex0_index, long vertex1_index)
    : VertexVertexCandidate(vertex0_index, vertex1_index)
{
}

VertexVertexFrictionConstraint::VertexVertexFrictionConstraint(
    const VertexVertexCandidate& candidate)
    : VertexVertexCandidate(candidate)
{
}

VertexVertexFrictionConstraint::VertexVertexFrictionConstraint(
    const VertexVertexConstraint& constraint)
    : VertexVertexCandidate(constraint.vertex0_index, constraint.vertex1_index)
    , multiplicity(constraint.multiplicity)
{
}

VectorMax12d VertexVertexFrictionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h) const
{
    VectorMax3d rel_u = relative_displacement(U);
    VectorMax2d tangent_rel_u =
        compute_potential_gradient_common(rel_u, epsv_times_h);
    tangent_rel_u *= multiplicity;
    return point_point_relative_mesh_displacements(
        tangent_rel_u, tangent_basis);
}

MatrixMax12d VertexVertexFrictionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h,
    const bool project_hessian_to_psd) const
{
    VectorMax3d rel_u = relative_displacement(U);
    auto TT = point_point_TT(tangent_basis);
    return compute_potential_hessian_common(
        rel_u, TT, epsv_times_h, project_hessian_to_psd, multiplicity);
}

///////////////////////////////////////////////////////////////////////////////

EdgeVertexFrictionConstraint::EdgeVertexFrictionConstraint(
    long edge_index, long vertex_index)
    : EdgeVertexCandidate(edge_index, vertex_index)
{
}

EdgeVertexFrictionConstraint::EdgeVertexFrictionConstraint(
    const EdgeVertexCandidate& candidate)
    : EdgeVertexCandidate(candidate)
{
}

EdgeVertexFrictionConstraint::EdgeVertexFrictionConstraint(
    const EdgeVertexConstraint& constraint)
    : EdgeVertexCandidate(constraint.edge_index, constraint.vertex_index)
    , multiplicity(constraint.multiplicity)
{
}

VectorMax12d EdgeVertexFrictionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h) const
{
    VectorMax3d rel_u = relative_displacement(U, E);
    VectorMax2d tangent_rel_u =
        compute_potential_gradient_common(rel_u, epsv_times_h);
    tangent_rel_u *= multiplicity;
    return point_edge_relative_mesh_displacements(
        tangent_rel_u, tangent_basis, closest_point[0]);
}

MatrixMax12d EdgeVertexFrictionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h,
    const bool project_hessian_to_psd) const
{
    VectorMax3d rel_u = relative_displacement(U, E);
    auto TT = point_edge_TT(tangent_basis, closest_point[0]);
    return compute_potential_hessian_common(
        rel_u, TT, epsv_times_h, project_hessian_to_psd, multiplicity);
}

///////////////////////////////////////////////////////////////////////////////

EdgeEdgeFrictionConstraint::EdgeEdgeFrictionConstraint(
    long edge0_index, long edge1_index)
    : EdgeEdgeCandidate(edge0_index, edge1_index)
{
}

EdgeEdgeFrictionConstraint::EdgeEdgeFrictionConstraint(
    const EdgeEdgeCandidate& candidate)
    : EdgeEdgeCandidate(candidate)
{
}

EdgeEdgeFrictionConstraint::EdgeEdgeFrictionConstraint(
    const EdgeEdgeConstraint& constraint)
    : EdgeEdgeCandidate(constraint.edge0_index, constraint.edge1_index)
{
}

VectorMax12d EdgeEdgeFrictionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h) const
{
    VectorMax3d rel_u = relative_displacement(U, E);
    VectorMax2d tangent_rel_u =
        compute_potential_gradient_common(rel_u, epsv_times_h);
    return edge_edge_relative_mesh_displacements(
        tangent_rel_u, tangent_basis, closest_point);
}

MatrixMax12d EdgeEdgeFrictionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h,
    const bool project_hessian_to_psd) const
{
    VectorMax3d rel_u = relative_displacement(U, E);
    auto TT = edge_edge_TT(tangent_basis, closest_point);
    return compute_potential_hessian_common(
        rel_u, TT, epsv_times_h, project_hessian_to_psd);
}

///////////////////////////////////////////////////////////////////////////////

FaceVertexFrictionConstraint::FaceVertexFrictionConstraint(
    long face_index, long vertex_index)
    : FaceVertexCandidate(face_index, vertex_index)
{
}

FaceVertexFrictionConstraint::FaceVertexFrictionConstraint(
    const FaceVertexCandidate& candidate)
    : FaceVertexCandidate(candidate)
{
}

FaceVertexFrictionConstraint::FaceVertexFrictionConstraint(
    const FaceVertexConstraint& constraint)
    : FaceVertexCandidate(constraint.face_index, constraint.vertex_index)
{
}

VectorMax12d FaceVertexFrictionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h) const
{
    VectorMax3d rel_u = relative_displacement(U, F);
    VectorMax2d tangent_rel_u =
        compute_potential_gradient_common(rel_u, epsv_times_h);
    return point_triangle_relative_mesh_displacements(
        tangent_rel_u, tangent_basis, closest_point);
}

MatrixMax12d FaceVertexFrictionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h,
    const bool project_hessian_to_psd) const
{
    VectorMax3d rel_u = relative_displacement(U, F);
    assert(rel_u.size() == 3);
    auto TT = point_triangle_TT(tangent_basis, closest_point);
    return compute_potential_hessian_common(
        rel_u, TT, epsv_times_h, project_hessian_to_psd);
}

///////////////////////////////////////////////////////////////////////////////

size_t FrictionConstraints::size() const
{
    return vv_constraints.size() + ev_constraints.size() + ee_constraints.size()
        + fv_constraints.size();
}

size_t FrictionConstraints::num_constraints() const
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

bool FrictionConstraints::empty() const
{
    return vv_constraints.empty() && ev_constraints.empty()
        && ee_constraints.empty() && fv_constraints.empty();
}

void FrictionConstraints::clear()
{
    vv_constraints.clear();
    ev_constraints.clear();
    ee_constraints.clear();
    fv_constraints.clear();
}

FrictionConstraint& FrictionConstraints::operator[](size_t idx)
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
    throw std::out_of_range("Friction constraint index is out of range!");
}

const FrictionConstraint& FrictionConstraints::operator[](size_t idx) const
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
    throw std::out_of_range("Friction constraint index is out of range!");
}

} // namespace ipc
