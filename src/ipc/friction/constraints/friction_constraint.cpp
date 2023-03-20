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
    const double epsv_times_h) const
{
    // μ N(xᵗ) f₀(‖u‖) (where u = T(xᵗ)ᵀv)

    // Compute u = PᵀΓv
    const VectorMax2d u = tangent_basis.transpose()
        * relative_velocity(dof(velocities, edges, faces));

    return weight * mu * normal_force_magnitude * f0_SF(u.norm(), epsv_times_h);
}

VectorMax12d FrictionConstraint::compute_potential_gradient(
    const Eigen::MatrixXd& velocities,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    double epsv_times_h) const
{
    assert(epsv_times_h > 0);
    // ∇ₓ μ N(xᵗ) f₀(‖u‖) (where u = T(xᵗ)ᵀv)
    //  = μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u

    // Compute u = PᵀΓv
    const VectorMax2d u = tangent_basis.transpose()
        * relative_velocity(dof(velocities, edges, faces));

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        relative_velocity_matrix().transpose() * tangent_basis;

    // Compute f₁(‖u‖)/‖u‖
    const double f1_over_norm_u = f1_SF_over_x(u.norm(), epsv_times_h);

    // μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u ∈ ℝⁿ
    // (n×2)(2×1) = (n×1)
    return T * ((weight * mu * normal_force_magnitude * f1_over_norm_u) * u);
}

MatrixMax12d FrictionConstraint::compute_potential_hessian(
    const Eigen::MatrixXd& velocities,
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
    const double epsv_times_h,
    bool project_hessian_to_psd) const
{
    assert(epsv_times_h > 0);
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
    const double f1_over_norm_u = f1_SF_over_x(norm_u, epsv_times_h);

    // Compute μ N(xᵗ)
    const double scale = weight * mu * normal_force_magnitude;

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
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
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

    const VectorMax12d x = dof(X, edges, faces);
    const VectorMax12d ut = dof(Ut, edges, faces);
    const VectorMax12d u = dof(U, edges, faces);

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
    const MatrixMax<double, 3, 12> Gamma = relative_velocity_matrix(beta);

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
    const Eigen::MatrixXi& edges,
    const Eigen::MatrixXi& faces,
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
    // Compute F
    assert(X.size() == U.size() && Ut.size() == U.size());
    int dim = U.cols();
    int n = dim * num_vertices();

    const VectorMax12d x = dof(X, edges, faces);
    const VectorMax12d ut = dof(Ut, edges, faces);
    const VectorMax12d u = dof(U, edges, faces);

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
    const MatrixMax<double, 3, 12> Gamma = relative_velocity_matrix(beta);

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
                relative_velocity_matrix_jacobian(beta);

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
    // Compute F = ∇(-μ N f₁(‖τ‖)/‖τ‖ T τ)
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
