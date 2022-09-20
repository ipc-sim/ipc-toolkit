#include "friction_constraint.hpp"

#include <ipc/friction/closest_point.hpp>
#include <ipc/friction/tangent_basis.hpp>
#include <ipc/friction/normal_force_magnitude.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_edge.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/point_triangle.hpp>

#include <ipc/config.hpp>

#include <stdexcept> // std::out_of_range

namespace ipc {

void FrictionConstraint::init(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double dhat,
    const double barrier_stiffness,
    const double dmin)
{
    tangent_basis.resize(V.cols(), V.cols() - 1); // do this to initialize dim()
    VectorMax12d x = select_dofs(V, E, F);
    closest_point = compute_closest_point(x);
    tangent_basis = compute_tangent_basis(x);
    normal_force_magnitude =
        compute_normal_force_magnitude(x, dhat, barrier_stiffness, dmin);
}

VectorMax12d FrictionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    double epsv_times_h) const
{
    assert(epsv_times_h > 0);
    // ∇ₓ μ N(xᵗ) f₀(‖u‖) (where u = T(xᵗ)ᵀ(x - xᵗ))
    //  = μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u

    // Compute u = PᵀΓ(x - xᵗ)
    const VectorMax2d u =
        tangent_basis.transpose() * relative_displacement(select_dofs(U, E, F));

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        relative_displacement_matrix().transpose() * tangent_basis;

    // Compute f₁(‖ū‖)/‖ū‖
    const double f1_over_norm_u = f1_SF_over_x(u.norm(), epsv_times_h);

    // μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u ∈ ℝⁿ
    return weight * mu * normal_force_magnitude * f1_over_norm_u
        * (T * u); // (n×2)(2×1) = (n×1)
}

MatrixMax12d FrictionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& U,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const double epsv_times_h,
    bool project_hessian_to_psd) const
{
    assert(epsv_times_h > 0);
    // ∇ₓ μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u (where u = T(xᵗ)ᵀ (x - xᵗ))
    //  = μ N T [(f₁'(‖u‖)‖u‖ − f₁(‖u‖))/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
    //  = μ N T [f₂(‖u‖) uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ

    // Compute u = PᵀΓ(x - xᵗ)
    const VectorMax2d u =
        tangent_basis.transpose() * relative_displacement(select_dofs(U, E, F));

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        relative_displacement_matrix().transpose() * tangent_basis;

    // Compute ‖u‖
    const double norm_u = u.norm();

    // Compute f₁(‖u‖)/‖u‖
    const double f1_over_norm_u = f1_SF_over_x(norm_u, epsv_times_h);

    // Compute μ N(xᵗ)
    double scale = weight * mu * normal_force_magnitude;

    MatrixMax12d hess;
    if (norm_u >= epsv_times_h) {
        // f₁(‖u‖) = 1 ⟹ f₁'(‖u‖) = 0
        //  ⟹ ∇²D(x) = μ N T [-f₁(‖u‖)/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
        //            = μ N T [-f₁(‖u‖)/‖u‖ uuᵀ/‖u‖² + f₁(‖u‖)/‖u‖ I] Tᵀ
        //            = μ N T [f₁(‖u‖)/‖u‖ (I - uuᵀ/‖u‖²)] Tᵀ
        //  ⟹ no SPD projection needed because f₁(‖u‖)/‖u‖ ≥ 0
        if (dim() == 2) {
            // I - uuᵀ/‖u‖² = 1 - u²/u² = 0 ⟹ ∇²D(x) = 0
            int n = T.rows(); // num vars
            hess.setZero(n, n);
        } else {
            assert(dim() == 3);
            // I - uuᵀ/‖u‖² = ūūᵀ / ‖u‖² (where ū⋅u = 0)
            Eigen::Vector2d u_perp(-u[1], u[0]);
            hess = // grouped to reduce number of operations
                (T * ((scale * f1_over_norm_u / (norm_u * norm_u)) * u_perp))
                * (u_perp.transpose() * T.transpose());
        }
    } else if (norm_u == 0) {
        // ∇²D = μ N T [(f₁'(‖u‖)‖u‖ − f₁(‖u‖))/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
        // lim_{‖u‖→0} ∇²D = μ N T [f₁(‖u‖)/‖u‖ I] Tᵀ
        // no SPD projection needed because μ N f₁(‖ū‖)/‖ū‖ ≥ 0
        hess = scale * f1_over_norm_u * T * T.transpose();
    } else {
        // ∇²D(x) = μ N T [f₂(‖u‖) uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
        //  ⟹ only need to project the inner 2x2 matrix to SPD
        const double f2 = df1_x_minus_f1_over_x3(norm_u, epsv_times_h);

        MatrixMax2d inner_hess = f2 * u * u.transpose();
        inner_hess.diagonal().array() += f1_over_norm_u;
        inner_hess *= scale;
        if (project_hessian_to_psd) {
            inner_hess = project_to_psd(inner_hess);
        }

        hess = T * inner_hess * T.transpose();
    }

    return hess;
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
    const double dmin,
    const bool no_mu) const
{
    // x is the rest position
    // uᵗ is the displacment at the begining of the timestep
    // uᵢ is the displacment at the begginging of the lagged solve
    // u is the displacment at the end of the timestep
    // (uᵢ = uᵗ; uₙ = u)
    //
    // Static simulation:
    // τ = T(x + uᵢ)ᵀu is the tangential displacment
    // F(x, u) = -μ N(x + uᵢ) f₁(‖τ‖)/‖τ‖ T(x + uᵢ) τ
    //
    // Time-dependent simulation:
    // τ = T(x + uᵢ)ᵀ(u - uᵗ) is the tangential displacment
    // F(x, uᵗ, u) = -μ N(x + uᵢ) f₁(‖τ‖)/‖τ‖ T(x + uᵢ) τ
    assert(X.size() == U.size() && Ut.size() == U.size());

    const VectorMax12d x = select_dofs(X, E, F);
    const VectorMax12d ut = select_dofs(Ut, E, F);
    const VectorMax12d u = select_dofs(U, E, F);

    // Assume uᵢ = uᵗ
    VectorMax12d x_plus_ui = x + ut;

    // Assume uᵢ = u
    // VectorMax12d x_plus_ui = x + u;

    // Compute N(x + uᵢ)
    double N = compute_normal_force_magnitude(
        x_plus_ui, dhat, barrier_stiffness, dmin);

    // Compute P
    const MatrixMax<double, 3, 2> P = compute_tangent_basis(x_plus_ui);

    // compute β
    const VectorMax2d beta = compute_closest_point(x_plus_ui);

    // Compute Γ
    const MatrixMax<double, 3, 12> Gamma = relative_displacement_matrix(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;

    // Compute τ = PᵀΓ(u - uᵗ)
    const VectorMax2d tau = T.transpose() * (u - ut);

    // Compute f₁(‖τ‖)/‖τ‖
    const double f1_over_norm_tau = f1_SF_over_x(tau.norm(), epsv_times_h);

    // F = -μ N f₁(‖τ‖)/‖τ‖ T τ
    // NOTE: no_mu -> leave mu out of this function (i.e., assuming mu = 1)
    return -weight * (no_mu ? 1.0 : mu) * N * f1_over_norm_tau * T * tau;
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
    // x is the rest position
    // uᵗ is the displacment at the begining of the timestep
    // uᵢ is the displacment at the begginging of the lagged solve
    // u is the displacment at the end of the timestep
    // (uᵢ = uᵗ; uₙ = u)
    //
    // Static simulation:
    // τ = T(x + uᵢ)ᵀu is the tangential displacment
    // F(x, u) = -μ N(x + uᵢ) f₁(‖τ‖)/‖τ‖ T(x + uᵢ) τ
    //
    // Time-dependent simulation:
    // τ = T(x + uᵢ)ᵀ(u - uᵗ) is the tangential displacment
    // F(x, uᵗ, u) = -μ N(x + uᵢ) f₁(‖τ‖)/‖τ‖ T(x + uᵢ) τ
    //
    // Compute ∇F
    assert(X.size() == U.size() && Ut.size() == U.size());
    int dim = U.cols();
    int n = dim * num_vertices();

    const VectorMax12d x = select_dofs(X, E, F);
    const VectorMax12d ut = select_dofs(Ut, E, F);
    const VectorMax12d u = select_dofs(U, E, F);

    // Assume uᵢ = uᵗ
    VectorMax12d x_plus_ui = x + ut;
    bool need_jac_N_or_T = wrt != DiffWRT::U;

    // Assume uᵢ = u
    // VectorMax12d x_plus_ui = x + u;
    // bool need_jac_N_or_T = wrt != DiffWRT::Ut;

    // Compute N
    double N = compute_normal_force_magnitude(
        x_plus_ui, dhat, barrier_stiffness, dmin);

    // Compute ∇N
    VectorMax12d grad_N;
    if (need_jac_N_or_T) {
        // ∇ₓN = ∇ᵤN
        grad_N = compute_normal_force_magnitude_gradient(
            x_plus_ui, dhat, barrier_stiffness, dmin);
        assert(grad_N.array().isFinite().all());
    }

    // Compute P
    const MatrixMax<double, 3, 2> P = compute_tangent_basis(x_plus_ui);

    // Compute β
    const VectorMax2d beta = compute_closest_point(x_plus_ui);

    // Compute Γ
    const MatrixMax<double, 3, 12> Gamma = relative_displacement_matrix(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;

    // Compute ∇T
    MatrixMax<double, 144, 2> jac_T;
    if (need_jac_N_or_T) {
        jac_T.resize(n * n, dim - 1);
        // ∇T = ∇(ΓᵀP) = ∇ΓᵀP + Γᵀ∇P
        const MatrixMax<double, 36, 2> jac_P =
            compute_tangent_basis_jacobian(x_plus_ui);
        for (int i = 0; i < n; i++) {
            // ∂T/∂xᵢ += Γᵀ ∂P/∂xᵢ
            jac_T.middleRows(i * n, n) =
                Gamma.transpose() * jac_P.middleRows(i * dim, dim);
        }

        // Vertex-vertex does not have a closest point
        if (beta.size()) {
            // ∇Γ(β) = ∇ᵦΓ∇β ∈ ℝ^{d×n×n} ≡ ℝ^{nd×n}
            const MatrixMax<double, 2, 12> jac_beta =
                compute_closest_point_jacobian(x_plus_ui);
            const MatrixMax<double, 6, 12> jac_Gamma_wrt_beta =
                relative_displacement_matrix_jacobian(beta);

            MatrixMax<double, 36, 12> jac_Gamma(n * dim, n);
            for (int k = 0; k < n; k++) {
                for (int i = 0; i < dim; i++) {
                    // ∂Γᵢⱼ/∂xₖ = ∂Γᵢⱼ/∂β ⋅ ∂β/∂xₖ
                    jac_Gamma.row(k * dim + i) =
                        jac_Gamma_wrt_beta(
                            Eigen::seqN(i, beta.size(), dim), Eigen::all)
                            .eval()
                            .transpose()
                        * jac_beta.col(k);
                }
            }

            for (int i = 0; i < n; i++) {
                // ∂T/∂xᵢ += ∂Γ/∂xᵢᵀ P
                jac_T.middleRows(i * n, n) +=
                    jac_Gamma.middleRows(i * dim, dim).transpose() * P;
            }
        }
    }

    // Compute τ = PᵀΓ(u - uᵗ)
    const VectorMax12d u_minus_ut = u - ut;
    const VectorMax2d tau = P.transpose() * Gamma * u_minus_ut;

    // Compute ∇τ = ∇T(x + uᵢ)ᵀ(u - uᵗ) + T(x + uᵢ)ᵀ∇(u - uᵗ)
    MatrixMax<double, 2, 12> jac_tau;
    if (need_jac_N_or_T) {
        jac_tau.resize(dim - 1, n);
        // Compute ∇T(x + uᵢ)ᵀ(u - uᵗ)
        for (int i = 0; i < n; i++) {
            jac_tau.col(i) =
                jac_T.middleRows(i * n, n).transpose() * (u_minus_ut);
        }
    } else {
        jac_tau.setZero(dim - 1, n);
    }
    switch (wrt) {
    case DiffWRT::X:
        break; // Tᵀ ∇ₓ(u - uᵗ) = 0
    case DiffWRT::Ut:
        jac_tau -= T.transpose(); // Tᵀ ∇_{uᵗ}(u - uᵗ) = -Tᵀ
        break;
    case DiffWRT::U:
        jac_tau += T.transpose(); // Tᵀ ∇ᵤ(u - uᵗ) = Tᵀ
    }

    // Compute f₁(‖τ‖)/‖τ‖
    const double tau_norm = tau.norm();
    const double f1_over_norm_tau = f1_SF_over_x(tau_norm, epsv_times_h);

    // Compute ∇(f₁(‖τ‖)/‖τ‖)
    VectorMax12d grad_f1_over_norm_tau;
    if (tau_norm == 0) {
        // lim_{x→0} f₂(x)x² = 0
        grad_f1_over_norm_tau.setZero(n);
    } else {
        // ∇ (f₁(‖τ‖)/‖τ‖) = (f₁'(‖τ‖)‖τ‖ - f₁(‖τ‖)) / ‖τ‖³ τᵀ ∇τ
        double f2 = df1_x_minus_f1_over_x3(tau_norm, epsv_times_h);
        assert(std::isfinite(f2));
        grad_f1_over_norm_tau = f2 * tau.transpose() * jac_tau;
    }

    // Premultiplied values
    const VectorMax12d T_times_tau = T * tau;

    ///////////////////////////////////////////////////////////////////////////
    // Compute ∇F = ∇(-μ N f₁(‖τ‖)/‖τ‖ T τ)
    MatrixMax12d J = MatrixMax12d::Zero(n, n);

    // = -μ f₁(‖τ‖)/‖τ‖ (T τ) [∇N]ᵀ
    if (need_jac_N_or_T) {
        J = f1_over_norm_tau * T_times_tau * grad_N.transpose();
    }

    // + -μ N T τ [∇(f₁(‖τ‖)/‖τ‖)]
    J += N * T_times_tau * grad_f1_over_norm_tau.transpose();

    // + -μ N f₁(‖τ‖)/‖τ‖ [∇T] τ
    if (need_jac_N_or_T) {
        const VectorMax2d scaled_tau = N * f1_over_norm_tau * tau;
        for (int i = 0; i < n; i++) {
            // ∂J/∂xᵢ = ∂T/∂xᵢ * τ
            J.col(i) += jac_T.middleRows(i * n, n) * scaled_tau;
        }
    }

    // + -μ N f₁(‖τ‖)/‖τ‖ T [∇τ]
    J += N * f1_over_norm_tau * T * jac_tau;

    // NOTE: ∇ₓw(x) is not local to the contact pair (i.e., it involves more
    // than the 4 contacting vertices), so we do not have enough information
    // here to compute the gradient. Instead this should be handled outside of
    // the function. For a simple multiplicitive model (∑ᵢ wᵢ Fᵢ) this can be
    // done easily.
    J *= -weight * mu;

    return J;
}

double FrictionConstraint::compute_normal_force_magnitude(
    const VectorMax12d& x,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    return ipc::compute_normal_force_magnitude(
        compute_distance(x), dhat, barrier_stiffness, dmin);
}

VectorMax12d FrictionConstraint::compute_normal_force_magnitude_gradient(
    const VectorMax12d& x,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    return ipc::compute_normal_force_magnitude_gradient(
        compute_distance(x), compute_distance_gradient(x), dhat,
        barrier_stiffness, dmin);
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
{
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

double
VertexVertexFrictionConstraint::compute_distance(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_point_distance(x.head(dim()), x.tail(dim()));
}

VectorMax12d VertexVertexFrictionConstraint::compute_distance_gradient(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax6d grad_d;
    point_point_distance_gradient(x.head(dim()), x.tail(dim()), grad_d);
    return grad_d;
}

MatrixMax<double, 3, 2> VertexVertexFrictionConstraint::compute_tangent_basis(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_point_tangent_basis(x.head(dim()), x.tail(dim()));
}

MatrixMax<double, 36, 2>
VertexVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_point_tangent_basis_jacobian(x.head(dim()), x.tail(dim()));
}

VectorMax2d VertexVertexFrictionConstraint::compute_closest_point(
    const VectorMax12d& x) const
{
    return VectorMax2d();
}

MatrixMax<double, 2, 12>
VertexVertexFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& x) const
{
    return MatrixMax<double, 2, 12>();
}

MatrixMax<double, 3, 12>
VertexVertexFrictionConstraint::relative_displacement_matrix(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 0);
    return point_point_relative_displacement_matrix(dim());
}

MatrixMax<double, 6, 12>
VertexVertexFrictionConstraint::relative_displacement_matrix_jacobian(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 0);
    return point_point_relative_displacement_matrix_jacobian(dim());
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
{
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

double
EdgeVertexFrictionConstraint::compute_distance(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_edge_distance(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()),
        PointEdgeDistanceType::P_E);
}

VectorMax12d EdgeVertexFrictionConstraint::compute_distance_gradient(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax9d grad_d;
    point_edge_distance_gradient(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()),
        PointEdgeDistanceType::P_E, grad_d);
    return grad_d;
}

MatrixMax<double, 3, 2>
EdgeVertexFrictionConstraint::compute_tangent_basis(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_edge_tangent_basis(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()));
}

MatrixMax<double, 36, 2>
EdgeVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_edge_tangent_basis_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()));
}

VectorMax2d
EdgeVertexFrictionConstraint::compute_closest_point(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax2d closest_point(1);
    closest_point[0] = point_edge_closest_point(
        x.head(dim()), x.segment(dim(), dim()), x.tail(dim()));
    return closest_point;
}

MatrixMax<double, 2, 12>
EdgeVertexFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_edge_closest_point_jacobian(
               x.head(dim()), x.segment(dim(), dim()), x.tail(dim()))
        .transpose();
}

MatrixMax<double, 3, 12>
EdgeVertexFrictionConstraint::relative_displacement_matrix(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 1);
    return point_edge_relative_displacement_matrix(dim(), closest_point[0]);
}

MatrixMax<double, 6, 12>
EdgeVertexFrictionConstraint::relative_displacement_matrix_jacobian(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 1);
    return point_edge_relative_displacement_matrix_jacobian(
        dim(), closest_point[0]);
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
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

double EdgeEdgeFrictionConstraint::compute_distance(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    // The distance type is known because mollified PP and PE were skipped.
    return edge_edge_distance(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()), EdgeEdgeDistanceType::EA_EB);
}

VectorMax12d EdgeEdgeFrictionConstraint::compute_distance_gradient(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax12d grad_d;
    // The distance type is known because mollified PP and PE were skipped.
    edge_edge_distance_gradient(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()), EdgeEdgeDistanceType::EA_EB, grad_d);
    return grad_d;
}

MatrixMax<double, 3, 2>
EdgeEdgeFrictionConstraint::compute_tangent_basis(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return edge_edge_tangent_basis(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 36, 2>
EdgeEdgeFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return edge_edge_tangent_basis_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

VectorMax2d
EdgeEdgeFrictionConstraint::compute_closest_point(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return edge_edge_closest_point(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 2, 12>
EdgeEdgeFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return edge_edge_closest_point_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 3, 12>
EdgeEdgeFrictionConstraint::relative_displacement_matrix(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 2);
    return edge_edge_relative_displacement_matrix(dim(), closest_point);
}

MatrixMax<double, 6, 12>
EdgeEdgeFrictionConstraint::relative_displacement_matrix_jacobian(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 2);
    return edge_edge_relative_displacement_matrix_jacobian(
        dim(), closest_point);
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
    this->weight = constraint.weight;
    this->weight_gradient = constraint.weight_gradient;
}

double
FaceVertexFrictionConstraint::compute_distance(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_distance(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()), PointTriangleDistanceType::P_T);
}

VectorMax12d FaceVertexFrictionConstraint::compute_distance_gradient(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    VectorMax12d grad_d;
    point_triangle_distance_gradient(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()), PointTriangleDistanceType::P_T, grad_d);
    return grad_d;
}

MatrixMax<double, 3, 2>
FaceVertexFrictionConstraint::compute_tangent_basis(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_tangent_basis(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 36, 2>
FaceVertexFrictionConstraint::compute_tangent_basis_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_tangent_basis_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

VectorMax2d
FaceVertexFrictionConstraint::compute_closest_point(const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_closest_point(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 2, 12>
FaceVertexFrictionConstraint::compute_closest_point_jacobian(
    const VectorMax12d& x) const
{
    assert(x.size() == ndof());
    return point_triangle_closest_point_jacobian(
        x.head(dim()), x.segment(dim(), dim()), x.segment(2 * dim(), dim()),
        x.tail(dim()));
}

MatrixMax<double, 3, 12>
FaceVertexFrictionConstraint::relative_displacement_matrix(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 2);
    return point_triangle_relative_displacement_matrix(dim(), closest_point);
}

MatrixMax<double, 6, 12>
FaceVertexFrictionConstraint::relative_displacement_matrix_jacobian(
    const VectorMax2d& closest_point) const
{
    assert(closest_point.size() == 2);
    return point_triangle_relative_displacement_matrix_jacobian(
        dim(), closest_point);
}

///////////////////////////////////////////////////////////////////////////////

size_t FrictionConstraints::size() const
{
    return vv_constraints.size() + ev_constraints.size() + ee_constraints.size()
        + fv_constraints.size();
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
