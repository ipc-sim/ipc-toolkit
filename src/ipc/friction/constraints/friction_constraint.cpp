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
    const Eigen::MatrixXd& positions,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const double barrier_stiffness,
    const double dmin)
{
    // do this to initialize dim()
    const int dim = positions.cols();
    tangent_basis.resize(dim, dim - 1);

    const VectorMax12d pos = dof(positions, edges, faces);
    closest_point = compute_closest_point(pos);
    tangent_basis = compute_tangent_basis(pos);
    normal_force_magnitude =
        compute_normal_force_magnitude(pos, dhat, barrier_stiffness, dmin);
}

double FrictionConstraint::compute_potential(
    const Eigen::MatrixXd& velocities,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double epsv) const
{
    // μ N(xᵗ) f₀(‖u‖) (where u = T(xᵗ)ᵀv)

    // Compute u = PᵀΓv
    const VectorMax2d u = tangent_basis.transpose()
        * relative_velocity(dof(velocities, edges, faces));

    return weight * mu * normal_force_magnitude * f0_SF(u.norm(), epsv);
}

VectorMax12d FrictionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& velocities,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double epsv) const
{
    assert(epsv > 0);
    // ∇ₓ μ N(xᵗ) f₀(‖u‖) (where u = T(xᵗ)ᵀv)
    //  = μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u

    // Compute u = PᵀΓv
    const VectorMax2d u = tangent_basis.transpose()
        * relative_velocity(dof(velocities, edges, faces));

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        relative_velocity_matrix().transpose() * tangent_basis;

    // Compute f₁(‖u‖)/‖u‖
    const double f1_over_norm_u = f1_SF_over_x(u.norm(), epsv);

    // μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u ∈ ℝⁿ
    // (n×2)(2×1) = (n×1)
    return T * ((weight * mu * normal_force_magnitude * f1_over_norm_u) * u);
}

MatrixMax12d FrictionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& velocities,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double epsv,
    bool project_hessian_to_psd) const
{
    assert(epsv > 0);
    // ∇ₓ μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u (where u = T(xᵗ)ᵀ v)
    //  = μ N T [(f₁'(‖u‖)‖u‖ − f₁(‖u‖))/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
    //  = μ N T [f₂(‖u‖) uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ

    // Compute u = PᵀΓv
    const VectorMax2d u = tangent_basis.transpose()
        * relative_velocity(dof(velocities, edges, faces));

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        relative_velocity_matrix().transpose() * tangent_basis;

    // Compute ‖u‖
    const double norm_u = u.norm();

    // Compute f₁(‖u‖)/‖u‖
    const double f1_over_norm_u = f1_SF_over_x(norm_u, epsv);

    // Compute μ N(xᵗ)
    const double scale = weight * mu * normal_force_magnitude;

    MatrixMax12d hess;
    if (norm_u >= epsv) {
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
            const Eigen::Vector2d u_perp(-u[1], u[0]);
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
        const double f2 = df1_x_minus_f1_over_x3(norm_u, epsv);

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
    const Eigen::MatrixXd& rest_positions,
    const Eigen::MatrixXd& lagged_displacements,
    const Eigen::MatrixXd& velocities,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const double barrier_stiffness,
    const double epsv,
    const double dmin,
    const bool no_mu) const
{
    // x is the rest position
    // u is the displacment at the begginging of the lagged solve
    // v is the current velocity
    //
    // τ = T(x + u)ᵀv is the tangential sliding velocity
    // F(x, u, v) = -μ N(x + u) f₁(‖τ‖)/‖τ‖ T(x + u) τ
    assert(rest_positions.size() == lagged_displacements.size());
    assert(lagged_displacements.size() == velocities.size());

    const VectorMax12d x = dof(rest_positions, edges, faces);
    const VectorMax12d u = dof(lagged_displacements, edges, faces);
    const VectorMax12d v = dof(velocities, edges, faces);
    const VectorMax12d x_plus_u = x + u;

    // Compute N(x + u)
    const double N =
        compute_normal_force_magnitude(x_plus_u, dhat, barrier_stiffness, dmin);

    // Compute P
    const MatrixMax<double, 3, 2> P = compute_tangent_basis(x_plus_u);

    // compute β
    const VectorMax2d beta = compute_closest_point(x_plus_u);

    // Compute Γ
    const MatrixMax<double, 3, 12> Gamma = relative_velocity_matrix(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;

    // Compute τ = PᵀΓv
    const VectorMax2d tau = T.transpose() * v;

    // Compute f₁(‖τ‖)/‖τ‖
    const double f1_over_norm_tau = f1_SF_over_x(tau.norm(), epsv);

    // F = -μ N f₁(‖τ‖)/‖τ‖ T τ
    // NOTE: no_mu -> leave mu out of this function (i.e., assuming mu = 1)
    return -weight * (no_mu ? 1.0 : mu) * N * f1_over_norm_tau * T * tau;
}

MatrixMax12d FrictionConstraint::compute_force_jacobian(
    const Eigen::MatrixXd& rest_positions,
    const Eigen::MatrixXd& lagged_displacements,
    const Eigen::MatrixXd& velocities,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double dhat,
    const double barrier_stiffness,
    const double epsv,
    const DiffWRT wrt,
    const double dmin) const
{
    // x is the rest position
    // u is the displacment at the begginging of the lagged solve
    // v is the current velocity
    //
    // τ = T(x + u)ᵀv is the tangential sliding velocity
    // F(x, u, v) = -μ N(x + u) f₁(‖τ‖)/‖τ‖ T(x + u) τ
    //
    // Compute ∇F
    assert(rest_positions.size() == lagged_displacements.size());
    assert(lagged_displacements.size() == velocities.size());
    int dim = rest_positions.cols();
    int n = dim * num_vertices();

    const VectorMax12d x = dof(rest_positions, edges, faces);
    const VectorMax12d u = dof(lagged_displacements, edges, faces);
    const VectorMax12d v = dof(velocities, edges, faces);

    const VectorMax12d x_plus_u = x + u;
    const bool need_jac_N_or_T = wrt != DiffWRT::VELOCITIES;

    // Compute N
    const double N =
        compute_normal_force_magnitude(x_plus_u, dhat, barrier_stiffness, dmin);

    // Compute ∇N
    VectorMax12d grad_N;
    if (need_jac_N_or_T) {
        // ∇ₓN = ∇ᵤN
        grad_N = compute_normal_force_magnitude_gradient(
            x_plus_u, dhat, barrier_stiffness, dmin);
        assert(grad_N.array().isFinite().all());
    }

    // Compute P
    const MatrixMax<double, 3, 2> P = compute_tangent_basis(x_plus_u);

    // Compute β
    const VectorMax2d beta = compute_closest_point(x_plus_u);

    // Compute Γ
    const MatrixMax<double, 3, 12> Gamma = relative_velocity_matrix(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;

    // Compute ∇T
    MatrixMax<double, 144, 2> jac_T;
    if (need_jac_N_or_T) {
        jac_T.resize(n * n, dim - 1);
        // ∇T = ∇(ΓᵀP) = ∇ΓᵀP + Γᵀ∇P
        const MatrixMax<double, 36, 2> jac_P =
            compute_tangent_basis_jacobian(x_plus_u);
        for (int i = 0; i < n; i++) {
            // ∂T/∂xᵢ += Γᵀ ∂P/∂xᵢ
            jac_T.middleRows(i * n, n) =
                Gamma.transpose() * jac_P.middleRows(i * dim, dim);
        }

        // Vertex-vertex does not have a closest point
        if (beta.size()) {
            // ∇Γ(β) = ∇ᵦΓ∇β ∈ ℝ^{d×n×n} ≡ ℝ^{nd×n}
            const MatrixMax<double, 2, 12> jac_beta =
                compute_closest_point_jacobian(x_plus_u);
            const MatrixMax<double, 6, 12> jac_Gamma_wrt_beta =
                relative_velocity_matrix_jacobian(beta);

            for (int k = 0; k < n; k++) {
                for (int b = 0; b < beta.size(); b++) {
                    jac_T.middleRows(k * n, n) +=
                        jac_Gamma_wrt_beta.transpose().middleCols(b * dim, dim)
                        * (jac_beta(b, k) * P);
                }
            }
        }
    }

    // Compute τ = PᵀΓv
    const VectorMax2d tau = P.transpose() * Gamma * v;

    // Compute ∇τ = ∇T(x + u)ᵀv + T(x + u)ᵀ∇v
    MatrixMax<double, 2, 12> jac_tau;
    if (need_jac_N_or_T) {
        jac_tau.resize(dim - 1, n);
        // Compute ∇T(x + u)ᵀv
        for (int i = 0; i < n; i++) {
            jac_tau.col(i) = jac_T.middleRows(i * n, n).transpose() * v;
        }
    } else {
        jac_tau = T.transpose(); // Tᵀ ∇ᵥv = Tᵀ
    }

    // Compute f₁(‖τ‖)/‖τ‖
    const double tau_norm = tau.norm();
    const double f1_over_norm_tau = f1_SF_over_x(tau_norm, epsv);

    // Compute ∇(f₁(‖τ‖)/‖τ‖)
    VectorMax12d grad_f1_over_norm_tau;
    if (tau_norm == 0) {
        // lim_{x→0} f₂(x)x² = 0
        grad_f1_over_norm_tau.setZero(n);
    } else {
        // ∇ (f₁(‖τ‖)/‖τ‖) = (f₁'(‖τ‖)‖τ‖ - f₁(‖τ‖)) / ‖τ‖³ τᵀ ∇τ
        double f2 = df1_x_minus_f1_over_x3(tau_norm, epsv);
        assert(std::isfinite(f2));
        grad_f1_over_norm_tau = f2 * tau.transpose() * jac_tau;
    }

    // Premultiplied values
    const VectorMax12d T_times_tau = T * tau;

    // ------------------------------------------------------------------------
    // Compute J = ∇F = ∇(-μ N f₁(‖τ‖)/‖τ‖ T τ)
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
    const VectorMax12d& positions,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    return ipc::compute_normal_force_magnitude(
        compute_distance(positions), dhat, barrier_stiffness, dmin);
}

VectorMax12d FrictionConstraint::compute_normal_force_magnitude_gradient(
    const VectorMax12d& positions,
    const double dhat,
    const double barrier_stiffness,
    const double dmin) const
{
    return ipc::compute_normal_force_magnitude_gradient(
        compute_distance(positions), compute_distance_gradient(positions), dhat,
        barrier_stiffness, dmin);
}

} // namespace ipc
