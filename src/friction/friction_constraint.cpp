#include <ipc/friction/friction_constraint.hpp>

#include <stdexcept> // std::out_of_range

#include <ipc/friction/tangent_basis.hpp>
#include <ipc/friction/normal_force_magnitude.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

namespace ipc {

VectorMax12d FrictionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double epsv_times_h) const
{
    // ∇ₓ μ N(xᵗ) f₀(‖ū‖) (where ū = T(xᵗ)ᵀ u(xᵗ, x))
    //  = μ N(xᵗ) f₁(‖ū‖)/‖ū‖ ūᵀ T(xᵗ)ᵀ ∇ₓu(xᵗ, x)

    // compute u(xᵗ, x)
    const VectorMax3d global_rel_u = relative_displacement(U, E, F);

    // compute ∇ₓ u(xᵗ, x)
    const MatrixMax<double, 3, 12> jac_global_rel_u =
        relative_displacement_jacobian(U, E, F);

    // compute ū
    const VectorMax2d tangent_rel_u = tangent_basis.transpose() * global_rel_u;

    // compute f₁(‖ū‖)/‖ū‖
    const double f1_over_norm_tangent_u =
        f1_SF_over_x(tangent_rel_u.norm(), epsv_times_h);

    // μ N(xᵗ) f₁(‖ū‖)/‖ū‖ ūᵀ T(xᵗ)ᵀ ∇ₓu(xᵗ, x)
    return multiplicity() * mu * normal_force_magnitude * f1_over_norm_tangent_u
        * tangent_rel_u.transpose() * tangent_basis.transpose()
        * jac_global_rel_u;
}

MatrixMax12d FrictionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h,
    bool project_hessian_to_psd) const
{
    int dim = U.cols();
    // ∇ₓ μ N(xᵗ) f₁(‖ū‖)/‖ū‖ ūᵀ T(xᵗ)ᵀ ∇ₓu(xᵗ, x) (where ū = T(xᵗ)ᵀ u(xᵗ, x))
    //  = μ N ∇uᵀ T [(f₁'(‖ū‖)‖ū‖ − f₁(‖ū‖))/‖ū‖³ ūūᵀ + f₁(‖ū‖)/‖ū‖ I] Tᵀ ∇u
    //  = μ N ∇uᵀ T [f₂(‖ū‖) ūūᵀ + f₁(‖ū‖)/‖ū‖ I] Tᵀ ∇u

    // compute u(xᵗ, x)
    const VectorMax3d global_rel_u = relative_displacement(U, E, F);

    // compute ∇ₓ u(xᵗ, x)
    const MatrixMax<double, 3, 12> jac_global_rel_u =
        relative_displacement_jacobian(U, E, F);

    // compute ū
    const VectorMax2d tangent_rel_u = tangent_basis.transpose() * global_rel_u;

    // compute TJᵤ = T(xᵗ)ᵀ ∇ₓu(xᵗ, x)
    const MatrixMax<double, 2, 12> TJ_u =
        tangent_basis.transpose() * jac_global_rel_u;

    // compute ‖ū‖
    const double norm_u = tangent_rel_u.norm();

    // compute f₁(‖ū‖)/‖ū‖
    const double f1_over_norm_ubar = f1_SF_over_x(norm_u, epsv_times_h);

    // compute μ N(xᵗ)
    double scale = multiplicity() * mu * normal_force_magnitude;

    MatrixMax12d local_hess;
    if (norm_u >= epsv_times_h) {
        // f₁(‖ū‖) = 1
        //  ⟹ ∇²D(x) = μ N ∇uᵀ T [-f₁(‖ū‖)/‖ū‖³ ūūᵀ + f₁(‖ū‖)/‖ū‖ I] Tᵀ ∇u
        //            = μ N ∇uᵀ T [-f₁(‖ū‖)/‖ū‖ ūūᵀ/‖ū‖² + f₁(‖ū‖)/‖ū‖ I] Tᵀ ∇u
        //            = μ N ∇uᵀ T [f₁(‖ū‖)/‖ū‖ (I - ūūᵀ/‖ū‖²)] Tᵀ ∇u
        //  ⟹ no SPD projection needed because f₁(‖ū‖)/‖ū‖ ≥ 0
        // ∧ dim = 2 ⟹
        if (dim == 2) {
            // I - ūūᵀ/‖ū‖² = 1 - ū²/ū² = 0 ⟹ ∇²D(x) = 0
            int n = jac_global_rel_u.cols(); // num_vars
            local_hess.setZero(n, n);
        } else {
            assert(dim == 3);
            // I - ūūᵀ/‖ū‖² = ūᵖ(ūᵖ)ᵀ / ‖ū‖² (where ūᵖ⋅ū = 0)
            Eigen::Vector2d u_perp(-tangent_rel_u[1], tangent_rel_u[0]);
            local_hess = // grouped to reduce number of operations
                (TJ_u.transpose()
                 * ((scale * f1_over_norm_ubar / (norm_u * norm_u)) * u_perp))
                * (u_perp.transpose() * TJ_u);
        }
    } else if (norm_u == 0) {
        // ∇²D = μ N ∇uᵀT[(f₁'(‖ū‖)‖ū‖ − f₁(‖ū‖))/‖ū‖³ ūūᵀ + f₁(‖ū‖)/‖ū‖ I]Tᵀ∇u
        // lim_{‖ū‖→0} ∇²D = μ N ∇uᵀT [f₁(‖ū‖)/‖ū‖ I] Tᵀ ∇u
        // no SPD projection needed because μ N f₁(‖ū‖)/‖ū‖ ≥ 0
        local_hess = scale * f1_over_norm_ubar * TJ_u.transpose() * TJ_u;
    } else {
        // ∇²D(x) = μ N ∇uᵀ T [f₂(‖ū‖) ūūᵀ + f₁(‖ū‖)/‖ū‖ I] Tᵀ ∇u
        //  ⟹ only need to project the inner 2x2 matrix to SPD
        const double f2 = df1_x_minus_f1_over_x3(norm_u, epsv_times_h);

        MatrixMax2d inner_hess = f2 * tangent_rel_u * tangent_rel_u.transpose();
        inner_hess.diagonal().array() += f1_over_norm_ubar;
        inner_hess *= scale;
        if (project_hessian_to_psd) {
            inner_hess = project_to_psd(inner_hess);
        }

        local_hess = TJ_u.transpose() * inner_hess * TJ_u;
    }

    return local_hess;
}

VectorMax12d FrictionConstraint::compute_force(
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const double dmin) const
{
    // X is the rest position of all vertices
    // Uᵗ is the displacment at the begining of the timestep of all vertices
    // Uⁱ is the displacment at the begginging of the lagged solve
    // U is the displacment at the end of the timestep of all vertices
    // (U⁰ = Uᵗ; Uⁿ = U)
    //
    // Static simulation:
    // u(U) is the displacment of the point of contact
    // ū = T(X + U)ᵀu(U) is the tangential displacment of u
    // F(X, U) = -μ N(X + U) f₁(‖ū‖)/‖ū‖ ūᵀ T(X + U)ᵀ ∇u(U)
    //
    // Time-dependent simulation:
    // u(U - Uᵗ) is the displacment of the point of contact over the timestep
    // ū = T(X + Uⁱ)ᵀu(U - Uᵗ) is the tangential displacment of u
    // F(X, Uᵗ, U) = -μ N(X + Uⁱ) f₁(‖ū‖)/‖ū‖ ūᵀ T(X + Uⁱ)ᵀ ∇ᵤu(U - Uᵗ)
    const bool is_time_dependent = Ut.size() != 0;

    // Assume Uⁱ = U

    // Eigen::MatrixXd displaced_X = X + (is_time_dependent ? Ut : U);
    Eigen::MatrixXd displaced_X = X + U;

    double N = compute_normal_force_magnitude(
        displaced_X, E, F, dhat, barrier_stiffness, dmin);

    MatrixMax<double, 3, 2> T = compute_tangent_basis(displaced_X, E, F);

    // compute u
    const VectorMax3d u =
        relative_displacement(is_time_dependent ? (U - Ut) : U, E, F);

    // compute ∇ᵤu
    const MatrixMax<double, 3, 12> jac_u =
        relative_displacement_jacobian(is_time_dependent ? (U - Ut) : U, E, F);

    // compute ū
    const VectorMax2d ubar = T.transpose() * u;

    // compute f₁(‖ū‖)/‖ū‖
    const double f1_over_norm_ubar = f1_SF_over_x(ubar.norm(), epsv_times_h);

    // μ N(x) f₁(‖ū‖)/‖ū‖ ūᵀ T(x)ᵀ ∇ₓu(xᵗ, x)
    return -multiplicity() * mu * N * f1_over_norm_ubar * ubar.transpose()
        * T.transpose() * jac_u;
}

MatrixMax12d FrictionConstraint::compute_force_jacobian(
    const Eigen::MatrixXd& X,
    const Eigen::MatrixXd& Ut,
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double epsv_times_h,
    const DiffWRT wrt,
    const double dmin) const
{
    int dim = U.cols();
    int n = dim * vertex_indices(E, F).size();
    const bool is_time_dependent = Ut.size() != 0;
    assert(wrt != DiffWRT::Ut || is_time_dependent);
    assert(
        !is_time_dependent || (U.rows() == Ut.rows() && U.cols() == Ut.cols()));

    // X is the rest position of all vertices
    // Uᵗ is the displacment at the begining of the timestep of all vertices
    // Uⁱ is the displacment at the begginging of the lagged solve
    // U is the displacment at the end of the timestep of all vertices
    // (U⁰ = Uᵗ; Uⁿ = U)
    //
    // Static simulation:
    // u(U) is the displacment of the point of contact
    // ū = T(X + U)ᵀu(U) is the tangential displacment of u
    // F(X, U) = -μ N(X + U) f₁(‖ū‖)/‖ū‖ ūᵀ T(X + U)ᵀ ∇u(U)
    //
    // Time-dependent simulation:
    // u(U - Uᵗ) is the displacment of the point of contact over the timestep
    // ū = T(X + Uⁱ)ᵀu(U - Uᵗ) is the tangential displacment of u
    // F(X, Uᵗ, U) = -μ N(X + Uⁱ) f₁(‖ū‖)/‖ū‖ ūᵀ T(X + Uⁱ)ᵀ ∇ᵤu(U - Uᵗ)
    //
    // Compute ∇F

    // Assume Uⁱ = U

    // Boolean for if we need to compute the derivative of N and T.
    // bool need_jac_N_or_T = !is_time_dependent || wrt != DiffWRT::U;
    bool need_jac_N_or_T = !is_time_dependent || wrt != DiffWRT::Ut;

    // Eigen::MatrixXd displaced_X = X + (is_time_dependent ? Ut : U);
    Eigen::MatrixXd displaced_X = X + U;

    // Compute N
    double N = compute_normal_force_magnitude(
        displaced_X, E, F, dhat, barrier_stiffness, dmin);

    // Compute ∇N
    VectorMax12d grad_N;
    if (need_jac_N_or_T) {
        // ∇ₓN = ∇ᵤN
        grad_N = compute_normal_force_magnitude_gradient(
            displaced_X, E, F, dhat, barrier_stiffness, dmin);
    }

    // Compute T
    MatrixMax<double, 3, 2> T = compute_tangent_basis(displaced_X, E, F);

    // Compute ∇T
    MatrixMax<double, 3, 24> jac_T;
    if (need_jac_N_or_T) {
        // ∇ₓT = ∇ᵤT
        jac_T = compute_tangent_basis_jacobian(displaced_X, E, F);
    }

    // Compute u
    VectorMax3d u =
        relative_displacement(is_time_dependent ? (U - Ut) : U, E, F);

    // Compute ∇ᵤu
    const MatrixMax<double, 3, 12> jac_u_wrt_U =
        relative_displacement_jacobian(is_time_dependent ? (U - Ut) : U, E, F);

    // Compute ū
    VectorMax2d ubar = T.transpose() * u;

    // Compute f₁(‖ū‖)/‖ū‖
    double f1_over_norm_ubar = f1_SF_over_x(ubar.norm(), epsv_times_h);

    // Premultiplied values
    const MatrixMax<double, 2, 12> T_Ju = T.transpose() * jac_u_wrt_U;
    const RowVectorMax12d ubar_T_Ju = ubar.transpose() * T_Ju;

    ///////////////////////////////////////////////////////////////////////////
    // Compute ∇ū
    MatrixMax<double, 2, 12> jac_ubar;
    if (need_jac_N_or_T) {
        MatrixMax<double, 24, 1> JT_u = jac_T.transpose() * u;
        jac_ubar =
            Eigen::Map<MatrixMax<double, 2, 12>>(JT_u.data(), dim - 1, n);
    } else {
        jac_ubar.setZero(dim - 1, n);
    }
    switch (wrt) {
    case DiffWRT::X:
        // ∇ₓū = ∇ₓT(X + U)ᵀu(U) or ∇ₓū = ∇ₓT(X + Uⁱ)ᵀu(U - Uᵗ)
        break;
    case DiffWRT::Ut:
        // ∇_{Uᵗ} ū = ∇_{Uᵗ} T(X + Uⁱ)ᵀu(U - Uᵗ) + T(X + Uⁱ)ᵀ ∇_{Uᵗ}u(U - Uᵗ)
        jac_ubar -= T_Ju; // ∇_{Uᵗ} u = -∇_{U} u
        break;
    case DiffWRT::U:
        // ∇ᵤū = ∇ᵤT(X+U)ᵀu(U) + T(X+U)ᵀ∇ᵤu(U) or ∇ᵤū = T(X+Uⁱ)ᵀ∇ᵤu(U-Uᵗ)
        jac_ubar += T_Ju;
    }

    ///////////////////////////////////////////////////////////////////////////
    // compute ∇F
    MatrixMax12d J = MatrixMax12d::Zero(n, n);

    // = -μ [∇N] f₁(‖ū‖)/‖ū‖ ūᵀ Tᵀ ∇ᵤu
    if (need_jac_N_or_T) {
        for (int i = 0; i < n; i++) {
            J.col(i) = grad_N(i) * f1_over_norm_ubar * ubar_T_Ju;
        }
    }

    // + -μ N [(f₁'(‖ū‖)‖ū‖ - f₁(‖ū‖))/‖ū‖³ (ūᵀ ∇ū)ᵀ] ūᵀ Tᵀ ∇ᵤu
    double f2 = df1_x_minus_f1_over_x3(ubar.norm(), epsv_times_h);
    J += (N * f2 * ubar.transpose() * jac_ubar).transpose() * ubar_T_Ju;

    // + -μ N f₁(‖ū‖)/‖ū‖ (∇ū)ᵀ Tᵀ ∇ᵤu
    J += N * f1_over_norm_ubar * jac_ubar.transpose() * T_Ju;

    // + -μ N f₁(‖ū‖)/‖ū‖ ūᵀ ∇Tᵀ ∇ᵤu
    if (need_jac_N_or_T) {
        for (int i = 0; i < n; i++) {
            J.col(i) += N * f1_over_norm_ubar * ubar.transpose()
                * jac_T.middleCols((dim - 1) * i, (dim - 1)).transpose()
                * jac_u_wrt_U;
        }
    }

    // + -μ N f₁(‖ū‖)/‖ū‖ ūᵀ Tᵀ (∇²u = 0).

    J *= -multiplicity() * mu;

    return J;
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
    , m_multiplicity(constraint.multiplicity)
{
}

double VertexVertexFrictionConstraint::compute_normal_force_magnitude(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    const auto& p0 = V.row(vertex0_index);
    const auto& p1 = V.row(vertex1_index);
    return ipc::compute_normal_force_magnitude(
        point_point_distance(p0, p1), dhat, barrier_stiffness, dmin);
}

VectorMax12d
VertexVertexFrictionConstraint::compute_normal_force_magnitude_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    const auto& p0 = V.row(vertex0_index);
    const auto& p1 = V.row(vertex1_index);
    VectorMax6d grad_d;
    point_point_distance_gradient(p0, p1, grad_d);
    return ipc::compute_normal_force_magnitude_gradient(
        point_point_distance(p0, p1), grad_d, dhat, barrier_stiffness, dmin);
}

MatrixMax<double, 3, 2> VertexVertexFrictionConstraint::compute_tangent_basis(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_point_tangent_basis(
        V.row(vertex0_index), V.row(vertex1_index));
}

MatrixMax<double, 3, 24>
VertexVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_point_tangent_basis_jacobian(
        V.row(vertex0_index), V.row(vertex1_index));
}

MatrixMax<double, 3, 12>
VertexVertexFrictionConstraint::relative_displacement_jacobian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_point_relative_displacement_jacobian(
        U.row(vertex0_index), U.row(vertex1_index));
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
    , m_multiplicity(constraint.multiplicity)
{
}

double EdgeVertexFrictionConstraint::compute_normal_force_magnitude(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    const auto& p = V.row(vertex_index);
    const auto& e0 = V.row(E(edge_index, 0));
    const auto& e1 = V.row(E(edge_index, 1));
    return ipc::compute_normal_force_magnitude(
        point_edge_distance(p, e0, e1), dhat, barrier_stiffness, dmin);
}

VectorMax12d
EdgeVertexFrictionConstraint::compute_normal_force_magnitude_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    const auto& p = V.row(vertex_index);
    const auto& e0 = V.row(E(edge_index, 0));
    const auto& e1 = V.row(E(edge_index, 1));
    VectorMax9d grad_d;
    point_edge_distance_gradient(p, e0, e1, grad_d);
    return ipc::compute_normal_force_magnitude_gradient(
        point_edge_distance(p, e0, e1), grad_d, dhat, barrier_stiffness, dmin);
}

MatrixMax<double, 3, 2> EdgeVertexFrictionConstraint::compute_tangent_basis(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_edge_tangent_basis(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)));
}

MatrixMax<double, 3, 24>
EdgeVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_edge_tangent_basis_jacobian(
        V.row(vertex_index), V.row(E(edge_index, 0)), V.row(E(edge_index, 1)));
}

MatrixMax<double, 3, 12>
EdgeVertexFrictionConstraint::relative_displacement_jacobian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_edge_relative_displacement_jacobian(
        U.row(vertex_index), U.row(E(edge_index, 0)), U.row(E(edge_index, 1)),
        closest_point[0]);
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

double EdgeEdgeFrictionConstraint::compute_normal_force_magnitude(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    const auto& ea0 = V.row(E(edge0_index, 0));
    const auto& ea1 = V.row(E(edge0_index, 1));
    const auto& eb0 = V.row(E(edge1_index, 0));
    const auto& eb1 = V.row(E(edge1_index, 1));
    return ipc::compute_normal_force_magnitude(
        edge_edge_distance(ea0, ea1, eb0, eb1), dhat, barrier_stiffness, dmin);
}

VectorMax12d
EdgeEdgeFrictionConstraint::compute_normal_force_magnitude_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    const auto& ea0 = V.row(E(edge0_index, 0));
    const auto& ea1 = V.row(E(edge0_index, 1));
    const auto& eb0 = V.row(E(edge1_index, 0));
    const auto& eb1 = V.row(E(edge1_index, 1));
    VectorMax12d grad_d;
    edge_edge_distance_gradient(ea0, ea1, eb0, eb1, grad_d);
    return ipc::compute_normal_force_magnitude_gradient(
        edge_edge_distance(ea0, ea1, eb0, eb1), grad_d, dhat, barrier_stiffness,
        dmin);
}

MatrixMax<double, 3, 2> EdgeEdgeFrictionConstraint::compute_tangent_basis(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return edge_edge_tangent_basis(
        V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
        V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)));
}

MatrixMax<double, 3, 24>
EdgeEdgeFrictionConstraint::compute_tangent_basis_jacobian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return edge_edge_tangent_basis_jacobian(
        V.row(E(edge0_index, 0)), V.row(E(edge0_index, 1)),
        V.row(E(edge1_index, 0)), V.row(E(edge1_index, 1)));
}

MatrixMax<double, 3, 12>
EdgeEdgeFrictionConstraint::relative_displacement_jacobian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return edge_edge_relative_displacement_jacobian(
        U.row(E(edge0_index, 0)), U.row(E(edge0_index, 1)),
        U.row(E(edge1_index, 0)), U.row(E(edge1_index, 1)), closest_point);
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

double FaceVertexFrictionConstraint::compute_normal_force_magnitude(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    const auto& p = V.row(vertex_index);
    const auto& t0 = V.row(F(face_index, 0));
    const auto& t1 = V.row(F(face_index, 1));
    const auto& t2 = V.row(F(face_index, 2));
    return ipc::compute_normal_force_magnitude(
        point_triangle_distance(p, t0, t1, t2), dhat, barrier_stiffness, dmin);
}

VectorMax12d
FaceVertexFrictionConstraint::compute_normal_force_magnitude_gradient(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    const auto& p = V.row(vertex_index);
    const auto& t0 = V.row(F(face_index, 0));
    const auto& t1 = V.row(F(face_index, 1));
    const auto& t2 = V.row(F(face_index, 2));
    VectorMax12d grad_d;
    point_triangle_distance_gradient(p, t0, t1, t2, grad_d);
    return ipc::compute_normal_force_magnitude_gradient(
        point_triangle_distance(p, t0, t1, t2), grad_d, dhat, barrier_stiffness,
        dmin);
}

MatrixMax<double, 3, 2> FaceVertexFrictionConstraint::compute_tangent_basis(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_triangle_tangent_basis(
        V.row(vertex_index), V.row(F(face_index, 0)), V.row(F(face_index, 1)),
        V.row(F(face_index, 2)));
}

MatrixMax<double, 3, 24>
FaceVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_triangle_tangent_basis_jacobian(
        V.row(vertex_index), V.row(F(face_index, 0)), V.row(F(face_index, 1)),
        V.row(F(face_index, 2)));
}

MatrixMax<double, 3, 12>
FaceVertexFrictionConstraint::relative_displacement_jacobian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F) const
{
    return point_triangle_relative_displacement_jacobian(
        U.row(vertex_index), U.row(F(face_index, 0)), U.row(F(face_index, 1)),
        U.row(F(face_index, 2)), closest_point);
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
        num_constraints += vv_constraint.multiplicity();
    }
    for (const auto& ev_constraint : ev_constraints) {
        num_constraints += ev_constraint.multiplicity();
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
