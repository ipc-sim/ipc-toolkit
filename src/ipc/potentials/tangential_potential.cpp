#include "friction_potential.hpp"

#include <ipc/utils/local_to_global.hpp>

#include <tbb/combinable.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace ipc {

namespace {

    /// Scalar μ for dissipative potential / force. Semi-implicit lagged
    /// anisotropy uses cached effective coefficients for one evaluation.
    /// Always uses
    /// mu_s_effective_lagged / mu_k_effective_lagged (isotropic scalars are
    /// copied there by
    /// TangentialCollisions::reset_lagged_anisotropic_friction_coefficients
    /// after build; matchstick values by
    /// update_lagged_anisotropic_friction_coefficients).
    void friction_mu_for_evaluation(
        const TangentialCollision& collision,
        const bool no_mu,
        double& mu_s_out,
        double& mu_k_out)
    {
        if (no_mu) {
            mu_s_out = 1.0;
            mu_k_out = 1.0;
        } else {
            mu_s_out = collision.mu_s_effective_lagged;
            mu_k_out = collision.mu_k_effective_lagged;
        }
    }

} // namespace

// -- Cumulative methods -------------------------------------------------------

Eigen::VectorXd TangentialPotential::force(
    const TangentialCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
    Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
    Eigen::ConstRef<Eigen::MatrixXd> velocities,
    const NormalPotential& normal_potential,
    const double dmin,
    const bool no_mu) const
{
    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(velocities.size());
    }

    const int dim = velocities.cols();
    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::combinable<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(velocities.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            Eigen::VectorXd& global_force = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto& collision = collisions[i];

                const VectorMax12d local_force = force(
                    collision, collision.dof(rest_positions, edges, faces),
                    collision.dof(lagged_displacements, edges, faces),
                    collision.dof(velocities, edges, faces), //
                    normal_potential, dmin, no_mu);

                const std::array<index_t, 4> vis =
                    collision.vertex_ids(mesh.edges(), mesh.faces());

                local_gradient_to_global_gradient(
                    local_force, vis, dim, global_force);
            }
        });

    return storage.combine([](const Eigen::VectorXd& a,
                              const Eigen::VectorXd& b) { return a + b; });
}

Eigen::SparseMatrix<double> TangentialPotential::force_jacobian(
    const TangentialCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
    Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
    Eigen::ConstRef<Eigen::MatrixXd> velocities,
    const NormalPotential& normal_potential,
    const DiffWRT wrt,
    const double dmin) const
{
    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(
            velocities.size(), velocities.size());
    }

    const int dim = velocities.cols();
    Eigen::ConstRef<Eigen::MatrixXi> edges = mesh.edges();
    Eigen::ConstRef<Eigen::MatrixXi> faces = mesh.faces();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& jac_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const TangentialCollision& collision = collisions[i];

                const MatrixMax12d local_force_jacobian = force_jacobian(
                    collision, collision.dof(rest_positions, edges, faces),
                    collision.dof(lagged_displacements, edges, faces),
                    collision.dof(velocities, edges, faces), //
                    normal_potential, wrt, dmin);

                const std::array<index_t, 4> vis =
                    collision.vertex_ids(mesh.edges(), mesh.faces());

                local_hessian_to_global_triplets(
                    local_force_jacobian, vis, dim, jac_triplets,
                    mesh.num_vertices());
            }
        });

    Eigen::SparseMatrix<double> jacobian(velocities.size(), velocities.size());
    for (const auto& local_jac_triplets : storage) {
        Eigen::SparseMatrix<double> local_jacobian(
            velocities.size(), velocities.size());
        local_jacobian.setFromTriplets(
            local_jac_triplets.begin(), local_jac_triplets.end());
        jacobian += local_jacobian;
    }

    // if wrt == X then compute ∇ₓ w(x)
    if (wrt == DiffWRT::REST_POSITIONS) {
        for (int i = 0; i < collisions.size(); i++) {
            const TangentialCollision& collision = collisions[i];
            assert(collision.weight_gradient.size() == rest_positions.size());
            if (collision.weight_gradient.size() != rest_positions.size()) {
                throw std::runtime_error(
                    "Shape derivative is not computed for friction collision!");
            }

            VectorMax12d local_force = force(
                collision, collision.dof(rest_positions, edges, faces),
                collision.dof(lagged_displacements, edges, faces),
                collision.dof(velocities, edges, faces), //
                normal_potential, dmin);
            assert(collision.weight != 0);
            local_force /= collision.weight;

            Eigen::SparseVector<double> force(rest_positions.size());
            force.reserve(local_force.size());
            local_gradient_to_global_gradient(
                local_force, collision.vertex_ids(edges, faces), dim, force);

            jacobian += force * collision.weight_gradient.transpose();
        }
    }

    return jacobian;
}
// -- Single collision methods -------------------------------------------------

double TangentialPotential::operator()(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMax12d> velocities) const
{
    //
    // μₛ = μₖ:
    //   μ N(xᵗ) f₀(‖u‖) (where u = T(xᵗ)ᵀv)
    //
    // μₛ != μₖ:
    //   N(xᵗ) g₀(‖u‖) (where g₀(x) = ∫ μ(x) f₀(x) dx)
    //

    // Compute u = PᵀΓv
    const VectorMax2d u = collision.tangent_basis.transpose()
        * collision.relative_velocity(velocities);

    // Apply anisotropic scaling: u_aniso = mu_aniso ⊙ u
    // Handle both 2D tangent space (3D sim) and 1D tangent space (2D sim)
    const VectorMax2d u_aniso =
        collision.mu_aniso.head(u.size()).cwiseProduct(u);

    const int tangent_dim = u.size();
    double mu_s, mu_k;
    friction_mu_for_evaluation(collision, false, mu_s, mu_k);

    return collision.weight * collision.normal_force_magnitude
        * mu_f0(u_aniso.norm(), mu_s, mu_k);
}

VectorMax12d TangentialPotential::gradient(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMax12d> velocities) const
{
    //
    // μₛ = μₖ:
    //   ∇ᵥ μ N(xᵗ) f₀(‖u‖) = μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u
    //
    // μₛ != μₖ:
    //   ∇ᵥ N(xᵗ) g₀(‖u‖)/‖u‖ T(xᵗ) u = N(xᵗ) g₁(‖u‖)/‖u‖ T(xᵗ) u
    //

    // Compute u = PᵀΓv
    const VectorMax2d u = collision.tangent_basis.transpose()
        * collision.relative_velocity(velocities);

    // Apply anisotropic scaling: u_aniso = mu_aniso ⊙ u
    // Handle both 2D tangent space (3D sim) and 1D tangent space (2D sim)
    const int tangent_dim = u.size();
    const VectorMax2d u_aniso =
        collision.mu_aniso.head(tangent_dim).cwiseProduct(u);

    // Compute T = ΓᵀP (then scale columns by mu_aniso for ∂u_aniso/∂u)
    const MatrixMax<double, 12, 2> T =
        collision.relative_velocity_jacobian().transpose()
        * collision.tangent_basis;

    double mu_s, mu_k;
    friction_mu_for_evaluation(collision, false, mu_s, mu_k);

    // Apply anisotropic scaling to T: T_aniso = T * diag(mu_aniso)
    // This accounts for ∂u_aniso/∂u = diag(mu_aniso) in the chain rule
    MatrixMax<double, 12, 2> T_aniso = T;
    T_aniso.col(0) *= collision.mu_aniso[0];
    if (tangent_dim > 1) {
        T_aniso.col(1) *= collision.mu_aniso[1];
    }

    const double norm_u_aniso = u_aniso.norm();
    // Semi-implicit lagged anisotropy: no ∂μ_eff/∂u_aniso (μ fixed for this
    // evaluation).
    VectorMax2d grad_wrt_u_aniso =
        mu_f1_over_x(norm_u_aniso, mu_s, mu_k) * u_aniso;

    grad_wrt_u_aniso *= collision.weight * collision.normal_force_magnitude;

    return T_aniso * grad_wrt_u_aniso;
}

MatrixMax12d TangentialPotential::hessian(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMax12d> velocities,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    //
    // μₛ = μₖ:
    //   ∇ᵥ μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u
    //    = μ N T [(f₁'(‖u‖)‖u‖ − f₁(‖u‖))/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
    //    = μ N T [f₂(‖u‖) uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
    //
    // μₛ != μₖ:
    //   ∇ᵥ N(xᵗ) g₁(‖u‖)/‖u‖ T(xᵗ) u = N T [g₂(‖u‖) uuᵀ + g₁(‖u‖)/‖u‖ I] Tᵀ
    //

    // Compute u = PᵀΓv
    const VectorMax2d u = collision.tangent_basis.transpose()
        * collision.relative_velocity(velocities);

    // Apply anisotropic scaling: u_aniso = mu_aniso ⊙ u
    // Handle both 2D tangent space (3D sim) and 1D tangent space (2D sim)
    const VectorMax2d u_aniso =
        collision.mu_aniso.head(u.size()).cwiseProduct(u);

    // Get tangent space dimension (1 for 2D sim, 2 for 3D sim)
    const int tangent_dim = u.size();

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        collision.relative_velocity_jacobian().transpose()
        * collision.tangent_basis;

    // Apply anisotropic scaling to T: T_aniso = T * diag(mu_aniso)
    // This accounts for ∂u_aniso/∂u = diag(mu_aniso) in the chain rule
    MatrixMax<double, 12, 2> T_aniso = T;
    T_aniso.col(0) *= collision.mu_aniso[0];
    if (tangent_dim > 1) {
        T_aniso.col(1) *= collision.mu_aniso[1];
    }

    // Compute ‖u_aniso‖
    const double norm_u = u_aniso.norm();

    double mu_s, mu_k;
    friction_mu_for_evaluation(collision, false, mu_s, mu_k);

    // Compute μ(‖u_aniso‖) f₁(‖u_aniso‖)/‖u_aniso‖
    const double mu_f1_over_norm_u = mu_f1_over_x(norm_u, mu_s, mu_k);

    // Compute N(xᵗ)
    const double scale = collision.weight * collision.normal_force_magnitude;

    MatrixMax12d hess;
    if (is_dynamic(norm_u)) {
        // f₁(‖u‖) = 1 ⟹ f₂(‖u‖) = 0
        //  ⟹ ∇²D(v) = μ N T [-f₁(‖u‖)/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
        //            = μ N T [-f₁(‖u‖)/‖u‖ uuᵀ/‖u‖² + f₁(‖u‖)/‖u‖ I] Tᵀ
        //            = μ N T [f₁(‖u‖)/‖u‖ (I - uuᵀ/‖u‖²)] Tᵀ
        //  ⟹ no PSD projection needed because f₁(‖u‖)/‖u‖ ≥ 0
        if (project_hessian_to_psd != PSDProjectionMethod::NONE && scale <= 0) {
            hess.setZero(collision.ndof(), collision.ndof()); // -PSD = NSD ⟹ 0
        } else if (collision.dim() == 2) {
            // I - uuᵀ/‖u‖² = 1 - u²/u² = 0 ⟹ ∇²D(v) = 0
            hess.setZero(collision.ndof(), collision.ndof());
        } else {
            assert(collision.dim() == 3);
            // I - u_aniso u_anisoᵀ/‖u_aniso‖² = ūūᵀ / ‖u_aniso‖²
            // (where ū⋅u_aniso = 0)
            const Eigen::Vector2d u_perp(-u_aniso[1], u_aniso[0]);
            hess = // grouped to reduce number of operations
                (T_aniso
                 * ((scale * mu_f1_over_norm_u / (norm_u * norm_u)) * u_perp))
                * (u_perp.transpose() * T_aniso.transpose());
        }
    } else if (norm_u == 0) {
        // ∇²D = μ N T [(f₂(‖u‖)‖u‖ − f₁(‖u‖))/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
        // lim_{‖u‖→0} ∇²D = μ N T [f₁(‖u‖)/‖u‖ I] Tᵀ
        // no PSD projection needed because μ N f₁(‖ū‖)/‖ū‖ ≥ 0
        if (project_hessian_to_psd != PSDProjectionMethod::NONE && scale <= 0) {
            hess.setZero(collision.ndof(), collision.ndof()); // -PSD = NSD ⟹ 0
        } else {
            hess = scale * mu_f1_over_norm_u * T_aniso * T_aniso.transpose();
        }
    } else {
        // ∇²D(v) = μ N T [f₂(‖u_aniso‖) u_aniso u_anisoᵀ +
        // f₁(‖u_aniso‖)/‖u_aniso‖ I] Tᵀ
        //  ⟹ only need to project the inner 2x2 matrix to PSD
        const double f2 = mu_f2_x_minus_mu_f1_over_x3(norm_u, mu_s, mu_k);

        MatrixMax2d inner_hess = f2 * u_aniso * u_aniso.transpose();
        inner_hess.diagonal().array() += mu_f1_over_norm_u;
        inner_hess *= scale; // NOTE: negative scaling will be projected out
        inner_hess = project_to_psd(inner_hess, project_hessian_to_psd);

        hess = T_aniso * inner_hess * T_aniso.transpose();
    }

    return hess;
}

VectorMax12d TangentialPotential::force(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMax12d> rest_positions,       // = x
    Eigen::ConstRef<VectorMax12d> lagged_displacements, // = u
    Eigen::ConstRef<VectorMax12d> velocities,           // = v
    const NormalPotential& normal_potential,
    const double dmin,
    const bool no_mu) const
{
    // x is the rest position
    // u is the displacement at the beginning of the lagged solve
    // v is the current velocity
    //
    // τ = T(x + u)ᵀv is the tangential sliding velocity
    // τ_aniso = mu_aniso ⊙ τ is the anisotropically-scaled velocity
    // F(x, u, v) = -μ N(x + u) f₁(‖τ_aniso‖)/‖τ_aniso‖ T(x + u) τ_aniso
    //
    // Combined anisotropic friction model:
    //   1. mu_aniso velocity scaling: τ_aniso = diag(mu_aniso) · τ
    //   2. Scalar μ_s, μ_k in the formulas below are always
    //      mu_s_effective_lagged / mu_k_effective_lagged (isotropic scalars are
    //      copied there on build; matchstick values from
    //      TangentialCollisions::update_lagged_anisotropic_friction_coefficients).
    //      ‖τ_aniso‖ and T τ_aniso use the current velocity here.
    assert(rest_positions.size() == lagged_displacements.size());
    assert(rest_positions.size() == velocities.size());

    // const VectorMax12d x = dof(rest_positions, edges, faces);
    // const VectorMax12d u = dof(lagged_displacements, edges, faces);
    // const VectorMax12d v = dof(velocities, edges, faces);

    // x:
    const VectorMax12d lagged_positions = rest_positions + lagged_displacements;

    // Compute N(x + u)
    const double N = normal_potential.force_magnitude(
        collision.compute_distance(lagged_positions), dmin);

    // Compute P
    const MatrixMax<double, 3, 2> P =
        collision.compute_tangent_basis(lagged_positions);

    // compute β
    const VectorMax2d beta = collision.compute_closest_point(lagged_positions);

    // Compute Γ
    const MatrixMax<double, 3, 12> Gamma =
        collision.relative_velocity_jacobian(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;

    // Compute τ = PᵀΓv
    const VectorMax2d tau = T.transpose() * velocities;

    // Get tangent space dimension (1 for 2D sim, 2 for 3D sim)
    const int tangent_dim = tau.size();

    // Always apply mu_aniso velocity scaling first: tau_aniso = mu_aniso ⊙ tau
    // Handle both 2D tangent space (3D sim) and 1D tangent space (2D sim)
    const VectorMax2d tau_aniso =
        collision.mu_aniso.head(tangent_dim).cwiseProduct(tau);

    double mu_s, mu_k;
    friction_mu_for_evaluation(collision, no_mu, mu_s, mu_k);

    // Compute μ(‖τ_aniso‖) f₁(‖τ_aniso‖)/‖τ_aniso‖
    const double tau_aniso_norm = tau_aniso.norm();
    const double mu_f1_over_norm_tau = mu_f1_over_x(tau_aniso_norm, mu_s, mu_k);

    // F = -μ N f₁(‖τ_aniso‖)/‖τ_aniso‖ T τ_aniso
    // NOTE: no_mu -> leave mu out of this function (i.e., assuming mu = 1)
    // NOTE: Direction-dependent (matchstick) μ uses lagged scalars on the
    //       collision; refresh with TangentialCollisions::
    //       update_lagged_anisotropic_friction_coefficients.
    return -collision.weight * N * mu_f1_over_norm_tau * T * tau_aniso;
}

MatrixMax12d TangentialPotential::force_jacobian(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMax12d> rest_positions,       // = x
    Eigen::ConstRef<VectorMax12d> lagged_displacements, // = u
    Eigen::ConstRef<VectorMax12d> velocities,           // = v
    const NormalPotential& normal_potential,
    const DiffWRT wrt,
    const double dmin) const
{
    // x is the rest position
    // u is the displacement at the beginning of the lagged solve
    // v is the current velocity
    //
    // τ = T(x + u)ᵀv is the tangential sliding velocity
    // τ_aniso = mu_aniso ⊙ τ is the anisotropically-scaled velocity
    // F(x, u, v) = -μ N(x + u) f₁(‖τ_aniso‖)/‖τ_aniso‖ T(x + u) τ_aniso
    //
    // Compute ∇F with combined anisotropic friction model:
    //   - mu_aniso velocity scaling is ALWAYS applied via jac_tau_aniso
    //   - Directional (matchstick) μ is lagged on the collision; no ∂μ/∂τ in
    //     this Jacobian (see friction_mu_for_evaluation / update_lagged_…).
    assert(rest_positions.size() == lagged_displacements.size());
    assert(lagged_displacements.size() == velocities.size());
    const int ndof = rest_positions.size();
    const int dim = ndof / collision.num_vertices();
    assert(ndof % collision.num_vertices() == 0);

    // const VectorMax12d x = dof(rest_positions, edges, faces);
    // const VectorMax12d u = dof(lagged_displacements, edges, faces);
    // const VectorMax12d v = dof(velocities, edges, faces);

    // x + u:
    const VectorMax12d lagged_positions = rest_positions + lagged_displacements;
    const bool need_jac_N_or_T = wrt != DiffWRT::VELOCITIES;

    // Compute N
    const double N = normal_potential.force_magnitude(
        collision.compute_distance(lagged_positions), dmin);

    // Compute ∇N
    VectorMax12d grad_N;
    if (need_jac_N_or_T) {
        // ∇ₓN = ∇ᵤN
        grad_N = normal_potential.force_magnitude_gradient(
            collision.compute_distance(lagged_positions),
            collision.compute_distance_gradient(lagged_positions), dmin);
        assert(grad_N.array().isFinite().all());
    }

    // Compute P
    const MatrixMax<double, 3, 2> P =
        collision.compute_tangent_basis(lagged_positions);

    // Compute β
    const VectorMax2d beta = collision.compute_closest_point(lagged_positions);

    // Compute Γ
    const MatrixMax<double, 3, 12> Gamma =
        collision.relative_velocity_jacobian(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;

    // Compute ∇T
    MatrixMax<double, 24, 12> jac_T;
    if (need_jac_N_or_T) {
        // ∇T = ∇(ΓᵀP) = ∇ΓᵀP + Γᵀ∇P

        // Compute Γᵀ∇P
        const MatrixMax<double, 6, 12> jac_P =
            collision.compute_tangent_basis_jacobian(lagged_positions);
        jac_T = (Gamma.transpose() * jac_P.reshaped(P.rows(), P.cols() * ndof))
                    .reshaped(T.size(), ndof);

        // Vertex-vertex does not have a closest point
        if (beta.size() != 0) {
            // Compute: ∇ΓᵀP = (∇ᵦΓ ∇β)ᵀ P
            const MatrixMax<double, 2, 12> dbeta_dx = // ∇β
                collision.compute_closest_point_jacobian(lagged_positions);
            const MatrixMax<double, 36, 2> dGamma_dbeta = // ∇ᵦΓ
                collision.relative_velocity_dx_dbeta(beta);

            // 1. Precompute dT/dβ = [dΓ/dβ]ᵀ P
            MatrixMax<double, 24, 2> dT_dbeta(T.size(), beta.size());
            for (int b = 0; b < beta.size(); ++b) {
                dT_dbeta.col(b) =
                    (dGamma_dbeta.col(b).reshaped(dim, ndof).transpose() * P)
                        .reshaped();
            }

            // 2. Apply chain rule: dT/dβ * dβ/dx
            jac_T += dT_dbeta * dbeta_dx;
        }
    }

    // Compute τ = PᵀΓv
    const VectorMax2d tau = P.transpose() * Gamma * velocities;

    // Get tangent space dimension (1 for 2D sim, 2 for 3D sim)
    const int tangent_dim = tau.size();

    // Always apply mu_aniso velocity scaling first: tau_aniso = mu_aniso ⊙ tau
    // Handle both 2D tangent space (3D sim) and 1D tangent space (2D sim)
    const VectorMax2d tau_aniso =
        collision.mu_aniso.head(tangent_dim).cwiseProduct(tau);

    // Compute ∇τ = ∇T(x + u)ᵀv + T(x + u)ᵀ∇v
    MatrixMax<double, 2, 12> jac_tau;
    if (need_jac_N_or_T) {
        jac_tau.resize(dim - 1, ndof);
        // Compute ∇T(x + u)ᵀv
        jac_tau = (jac_T.reshaped(ndof, T.size()).transpose() * velocities)
                      .reshaped(jac_tau.rows(), jac_tau.cols());
    } else {
        jac_tau = T.transpose(); // Tᵀ ∇ᵥv = Tᵀ
    }

    // Compute ∇tau_aniso = diag(mu_aniso) * ∇tau (chain rule for mu_aniso
    // scaling)
    MatrixMax<double, 2, 12> jac_tau_aniso = jac_tau;
    jac_tau_aniso.row(0) *= collision.mu_aniso[0];
    if (tangent_dim > 1) {
        jac_tau_aniso.row(1) *= collision.mu_aniso[1];
    }

    double mu_s, mu_k;
    friction_mu_for_evaluation(collision, false, mu_s, mu_k);

    // Compute μ f₁(‖τ_aniso‖)/‖τ_aniso‖
    const double tau_aniso_norm = tau_aniso.norm();
    const double mu_f1_over_norm_tau = mu_f1_over_x(tau_aniso_norm, mu_s, mu_k);

    // ∇(μ f₁/‖τ‖) with μ lagged w.r.t. velocity (no ∂μ/∂τ in Jacobian w.r.t.
    // v).
    VectorMax12d grad_mu_f1_over_norm_tau;
    if (tau_aniso_norm == 0) {
        // lim_{x→0} f₂(x)x² = 0
        grad_mu_f1_over_norm_tau.setZero(ndof);
    } else {
        const double f2 =
            mu_f2_x_minus_mu_f1_over_x3(tau_aniso_norm, mu_s, mu_k);
        assert(std::isfinite(f2));
        grad_mu_f1_over_norm_tau = f2 * tau_aniso.transpose() * jac_tau_aniso;
    }

    // Premultiplied values
    const VectorMax12d T_times_tau = T * tau_aniso;

    // ------------------------------------------------------------------------
    // Compute J = ∇F = ∇(-μ N f₁(‖τ_aniso‖)/‖τ_aniso‖ T τ_aniso)
    MatrixMax12d J = MatrixMax12d::Zero(ndof, ndof);

    // = -μ f₁(‖τ_aniso‖)/‖τ_aniso‖ (T τ_aniso) [∇N]ᵀ
    if (need_jac_N_or_T) {
        J = mu_f1_over_norm_tau * T_times_tau * grad_N.transpose();
    }

    // + -N T τ_aniso [∇(μ f₁(‖τ_aniso‖)/‖τ_aniso‖)]
    J += N * T_times_tau * grad_mu_f1_over_norm_tau.transpose();

    // + -μ N f₁(‖τ_aniso‖)/‖τ_aniso‖ [∇T] τ_aniso
    if (need_jac_N_or_T) {
        const VectorMax2d scaled_tau = N * mu_f1_over_norm_tau * tau_aniso;
        for (int i = 0; i < ndof; i++) {
            // ∂J/∂xᵢ = ∂T/∂xᵢ * τ_aniso
            J.col(i) += jac_T.col(i).reshaped(T.rows(), T.cols()) * scaled_tau;
        }
    }

    // + -μ N f₁(‖τ_aniso‖)/‖τ_aniso‖ T [∇τ_aniso]
    J += N * mu_f1_over_norm_tau * T * jac_tau_aniso;

    // NOTE: ∇ₓw(x) is not local to the collision pair (i.e., it involves more
    // than the 4 colliding vertices), so we do not have enough information
    // here to compute the gradient. Instead this should be handled outside of
    // the function. For a simple multiplicative model (∑ᵢ wᵢ Fᵢ) this can be
    // done easily.
    J *= -collision.weight;

    return J;
}

Eigen::VectorXd TangentialPotential::smooth_contact_force(
    const TangentialCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
    Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
    Eigen::ConstRef<Eigen::MatrixXd> velocities,
    const double dmin,
    const bool no_mu) const
{
    if (collisions.empty()) {
        return Eigen::VectorXd::Zero(velocities.size());
    }

    const int dim = velocities.cols();
    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::enumerable_thread_specific<Eigen::VectorXd> storage(
        Eigen::VectorXd::Zero(velocities.size()));

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            Eigen::VectorXd& global_force = storage.local();
            for (size_t i = r.begin(); i < r.end(); i++) {
                const auto& collision = collisions[i];

                const VectorMaxNd local_force = smooth_contact_force(
                    collision, collision.dof(rest_positions, edges, faces),
                    collision.dof(lagged_displacements, edges, faces),
                    collision.dof(velocities, edges, faces), no_mu);

                const auto vis =
                    collision.vertex_ids(mesh.edges(), mesh.faces());

                local_gradient_to_global_gradient(
                    local_force, vis, dim, global_force);
            }
        });

    return storage.combine([](const Eigen::VectorXd& a,
                              const Eigen::VectorXd& b) { return a + b; });
}

Eigen::SparseMatrix<double> TangentialPotential::smooth_contact_force_jacobian(
    const TangentialCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
    Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
    Eigen::ConstRef<Eigen::MatrixXd> velocities,
    const SmoothContactParameters& params,
    const DiffWRT wrt,
    const double dmin,
    const bool no_mu) const
{
    if (collisions.empty()) {
        return Eigen::SparseMatrix<double>(
            velocities.size(), velocities.size());
    }

    const int dim = velocities.cols();
    const Eigen::MatrixXi& edges = mesh.edges();
    const Eigen::MatrixXi& faces = mesh.faces();

    tbb::enumerable_thread_specific<std::vector<Eigen::Triplet<double>>>
        storage;

    const Eigen::MatrixXd lagged_positions =
        rest_positions + lagged_displacements;

    tbb::parallel_for(
        tbb::blocked_range<size_t>(size_t(0), collisions.size()),
        [&](const tbb::blocked_range<size_t>& r) {
            auto& jac_triplets = storage.local();

            for (size_t i = r.begin(); i < r.end(); i++) {
                const TangentialCollision& collision = collisions[i];

                // This jacobian doesn't include the derivatives of normal
                // contact force
                const MatrixMaxNd local_force_jacobian =
                    collision.normal_force_magnitude
                    * smooth_contact_force_jacobian_unit(
                        collision,
                        collision.dof(lagged_positions, edges, faces),
                        collision.dof(velocities, edges, faces), wrt, false);

                const auto vis =
                    collision.vertex_ids(mesh.edges(), mesh.faces());

                local_hessian_to_global_triplets(
                    local_force_jacobian, vis, dim, jac_triplets,
                    mesh.num_vertices());

                if (wrt == DiffWRT::VELOCITIES) {
                    continue;
                }

                // The term that includes derivatives of normal contact force
                const VectorMaxNd local_force = smooth_contact_force(
                    collision, collision.dof(rest_positions, edges, faces),
                    collision.dof(lagged_displacements, edges, faces),
                    collision.dof(velocities, edges, faces), false, true);

                // normal_force_grad is the gradient of contact force norm
                Eigen::VectorXd normal_force_grad;
                std::vector<index_t> cc_vert_ids;
                Eigen::MatrixXd Xt = rest_positions + lagged_displacements;
                auto cc = collision.smooth_collision;
                const Eigen::VectorXd contact_grad =
                    cc->gradient(cc->dof(Xt), params);
                const Eigen::MatrixXd contact_hess =
                    cc->hessian(cc->dof(Xt), params);
                normal_force_grad =
                    (1 / contact_grad.norm()) * (contact_hess * contact_grad);
                cc_vert_ids = cc->vertex_ids();

                local_jacobian_to_global_triplets(
                    local_force * normal_force_grad.transpose(), vis,
                    cc_vert_ids, dim, jac_triplets, mesh.num_vertices(),
                    mesh.num_vertices());
            }
        });

    Eigen::SparseMatrix<double> jacobian(velocities.size(), velocities.size());
    for (const auto& local_jac_triplets : storage) {
        Eigen::SparseMatrix<double> local_jacobian(
            velocities.size(), velocities.size());
        local_jacobian.setFromTriplets(
            local_jac_triplets.begin(), local_jac_triplets.end());
        jacobian += local_jacobian;
    }

    return jacobian;
}

TangentialPotential::VectorMaxNd TangentialPotential::smooth_contact_force(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMaxNd> rest_positions,       // = x
    Eigen::ConstRef<VectorMaxNd> lagged_displacements, // = u
    Eigen::ConstRef<VectorMaxNd> velocities,           // = v
    const bool no_mu,
    const bool no_contact_force_multiplier) const
{
    // x is the rest position
    // u is the displacment at the begginging of the lagged solve
    // v is the current velocity
    //
    // τ = T(x + u)ᵀv is the tangential sliding velocity
    // F(x, u, v) = -μ N(x + u) f₁(‖τ‖)/‖τ‖ T(x + u) τ
    assert(rest_positions.size() == lagged_displacements.size());
    assert(rest_positions.size() == velocities.size());

    // const VectorMaxNd x = dof(rest_positions, edges, faces);
    // const VectorMax<double, STENCIL_NDOF> u =
    //     dof(lagged_displacements, edges, faces);
    // const VectorMaxNd v = dof(velocities, edges, faces);
    const VectorMaxNd lagged_positions =
        rest_positions + lagged_displacements; // = x

    // Compute N(x + u)
    // const double N = collision.compute_normal_force_magnitude(
    //     lagged_positions, barrier_potential, dmin);
    const double N = collision.normal_force_magnitude;

    // Compute P
    const MatrixMax<double, 3, 2> P =
        collision.compute_tangent_basis(lagged_positions);

    // compute β
    const VectorMax2d beta = collision.compute_closest_point(lagged_positions);

    // Compute Γ
    const MatrixMax<double, 3, STENCIL_NDOF> Gamma =
        collision.relative_velocity_jacobian(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, STENCIL_NDOF, 2> T = Gamma.transpose() * P;

    // Compute τ = PᵀΓv
    const VectorMax2d tau = T.transpose() * velocities;

    // Apply anisotropic scaling: tau_aniso = mu_aniso ⊙ tau
    // Handle both 2D tangent space (3D sim) and 1D tangent space (2D sim)
    const VectorMax2d tau_aniso =
        collision.mu_aniso.head(tau.size()).cwiseProduct(tau);

    // Get tangent space dimension (1 for 2D sim, 2 for 3D sim)
    const int tangent_dim = tau.size();

    double mu_s, mu_k;
    friction_mu_for_evaluation(collision, no_mu, mu_s, mu_k);

    const double mu_f1_over_norm_tau =
        mu_f1_over_x(tau_aniso.norm(), mu_s, mu_k);

    // F = -μ N f₁(‖tau_aniso‖)/‖tau_aniso‖ T tau_aniso
    // NOTE: no_mu -> leave mu out of this function (i.e., assuming mu = 1)
    return -collision.weight * (no_contact_force_multiplier ? 1.0 : N)
        * mu_f1_over_norm_tau * T * tau_aniso;
}

TangentialPotential::MatrixMaxNd
TangentialPotential::smooth_contact_force_jacobian_unit(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMaxNd> lagged_positions, // = x + u^t
    Eigen::ConstRef<VectorMaxNd> velocities,       // = v
    const DiffWRT wrt,
    const bool no_mu) const
{
    // x is the rest position
    // u is the displacment at the begginging of the lagged solve
    // v is the current velocity
    //
    // τ = T(x + u)ᵀv is the tangential sliding velocity
    // F(x, u, v) = -μ f₁(‖τ‖)/‖τ‖ T(x + u) τ
    //
    // Compute ∇F
    assert(lagged_positions.size() == velocities.size());
    const int ndof = lagged_positions.size();
    const int dim = ndof / collision.num_vertices();
    assert(ndof % collision.num_vertices() == 0);

    const bool need_jac_N_or_T = wrt != DiffWRT::VELOCITIES;

    // Compute P
    const MatrixMax<double, 3, 2> P =
        collision.compute_tangent_basis(lagged_positions);

    // Compute β
    const VectorMax2d beta = collision.compute_closest_point(lagged_positions);

    // Compute Γ
    const MatrixMax<double, 3, STENCIL_NDOF> Gamma =
        collision.relative_velocity_jacobian(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, STENCIL_NDOF, 2> T = Gamma.transpose() * P;

    // Compute ∇T
    MatrixMax<double, 2 * STENCIL_NDOF, STENCIL_NDOF> jac_T;
    if (need_jac_N_or_T) {
        // ∇T = ∇(ΓᵀP) = ∇ΓᵀP + Γᵀ∇P

        // Compute Γᵀ∇P
        const MatrixMax<double, 6, 12> jac_P =
            collision.compute_tangent_basis_jacobian(lagged_positions);
        jac_T = (Gamma.transpose() * jac_P.reshaped(P.rows(), P.cols() * ndof))
                    .reshaped(T.size(), ndof);

        // Vertex-vertex does not have a closest point
        if (beta.size()) {
            // Compute: ∇ΓᵀP = (∇ᵦΓ ∇β)ᵀ P
            const MatrixMax<double, 2, STENCIL_NDOF> dbeta_dx = // ∇β
                collision.compute_closest_point_jacobian(lagged_positions);
            const MatrixMax<double, 3 * STENCIL_NDOF, 2> dGamma_dbeta = // ∇ᵦΓ
                collision.relative_velocity_dx_dbeta(beta);

            // 1. Precompute dT/dβ = [dΓ/dβ]ᵀ P
            MatrixMax<double, 2 * STENCIL_NDOF, 2> dT_dbeta(
                T.size(), beta.size());
            for (int b = 0; b < beta.size(); ++b) {
                dT_dbeta.col(b) =
                    (dGamma_dbeta.col(b).reshaped(dim, ndof).transpose() * P)
                        .reshaped();
            }

            // 2. Apply chain rule: dT/dβ * dβ/dx
            jac_T += dT_dbeta * dbeta_dx;
        }
    }

    // Compute τ = PᵀΓv
    const VectorMax2d tau = P.transpose() * Gamma * velocities;

    // Get tangent space dimension (1 for 2D sim, 2 for 3D sim)
    const int tangent_dim = tau.size();

    // Apply anisotropic scaling: tau_aniso = mu_aniso ⊙ tau
    // Handle both 2D tangent space (3D sim) and 1D tangent space (2D sim)
    const VectorMax2d tau_aniso =
        collision.mu_aniso.head(tangent_dim).cwiseProduct(tau);

    // Compute ∇τ = ∇T(x + u)ᵀv + T(x + u)ᵀ∇v
    MatrixMax<double, 2, STENCIL_NDOF> jac_tau;
    if (need_jac_N_or_T) {
        jac_tau.resize(dim - 1, ndof);
        // Compute ∇T(x + u)ᵀv
        jac_tau = (jac_T.reshaped(ndof, T.size()).transpose() * velocities)
                      .reshaped(jac_tau.rows(), jac_tau.cols());
    } else {
        jac_tau = T.transpose(); // Tᵀ ∇ᵥv = Tᵀ
    }

    // Compute ∇tau_aniso = diag(mu_aniso) * ∇tau
    MatrixMax<double, 2, STENCIL_NDOF> jac_tau_aniso = jac_tau;
    jac_tau_aniso.row(0) *= collision.mu_aniso[0];
    if (tangent_dim > 1) {
        jac_tau_aniso.row(1) *= collision.mu_aniso[1];
    }

    const double tau_norm = tau_aniso.norm();
    double mu_s, mu_k;
    friction_mu_for_evaluation(collision, no_mu, mu_s, mu_k);
    const double f1_over_norm_tau = mu_f1_over_x(tau_norm, mu_s, mu_k);

    // Compute ∇(f₁(‖tau_aniso‖)/‖tau_aniso‖)
    VectorMaxNd grad_f1_over_norm_tau;
    if (tau_norm == 0) {
        // lim_{x→0} f₂(x)x² = 0
        grad_f1_over_norm_tau.setZero(ndof);
    } else {
        // ∇ (f₁(‖tau_aniso‖)/‖tau_aniso‖) = (f₁'(‖tau_aniso‖)‖tau_aniso‖ -
        // f₁(‖tau_aniso‖)) / ‖tau_aniso‖³ tau_anisoᵀ ∇tau_aniso
        double f2 = mu_f2_x_minus_mu_f1_over_x3(tau_norm, mu_s, mu_k);
        assert(std::isfinite(f2));
        grad_f1_over_norm_tau = f2 * tau_aniso.transpose() * jac_tau_aniso;
    }

    // Premultiplied values
    const VectorMaxNd T_times_tau = T * tau_aniso;

    // ------------------------------------------------------------------------
    // Compute J = ∇F = ∇(-μ N f₁(‖τ‖)/‖τ‖ T τ)
    MatrixMaxNd J = MatrixMaxNd::Zero(ndof, ndof);

    // + -μ N T τ [∇(f₁(‖τ‖)/‖τ‖)]
    J += T_times_tau * grad_f1_over_norm_tau.transpose();

    // + -μ N f₁(‖tau_aniso‖)/‖tau_aniso‖ [∇T] tau_aniso
    if (need_jac_N_or_T) {
        const VectorMax2d scaled_tau = f1_over_norm_tau * tau_aniso;
        for (int i = 0; i < ndof; i++) {
            // ∂J/∂xᵢ = ∂T/∂xᵢ * τ_aniso
            J.col(i) += jac_T.col(i).reshaped(T.rows(), T.cols()) * scaled_tau;
        }
    }

    // + -μ N f₁(‖tau_aniso‖)/‖tau_aniso‖ T [∇tau_aniso]
    J += f1_over_norm_tau * T * jac_tau_aniso;

    // NOTE: ∇ₓw(x) is not local to the collision pair (i.e., it involves more
    // than the 4 colliding vertices), so we do not have enough information
    // here to compute the gradient. Instead this should be handled outside of
    // the function. For a simple multiplicative model (∑ᵢ wᵢ Fᵢ) this can be
    // done easily.
    J *= -collision.weight;

    return J;
}

} // namespace ipc