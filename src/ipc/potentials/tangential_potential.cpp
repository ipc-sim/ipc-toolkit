#include "friction_potential.hpp"

#include <ipc/friction/smooth_mu.hpp>
#include <ipc/utils/local_to_global.hpp>

#include <tbb/combinable.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace ipc {

// -- Cumulative methods -------------------------------------------------------

Eigen::VectorXd TangentialPotential::force(
    const TangentialCollisions& collisions,
    const CollisionMesh& mesh,
    Eigen::ConstRef<Eigen::MatrixXd> rest_positions,
    Eigen::ConstRef<Eigen::MatrixXd> lagged_displacements,
    Eigen::ConstRef<Eigen::MatrixXd> velocities,
    const NormalPotential& normal_potential,
    const double normal_stiffness,
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
                    normal_potential, normal_stiffness, dmin, no_mu);

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
    const double normal_stiffness,
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
                    normal_potential, normal_stiffness, wrt, dmin);

                const std::array<index_t, 4> vis =
                    collision.vertex_ids(mesh.edges(), mesh.faces());

                local_hessian_to_global_triplets(
                    local_force_jacobian, vis, dim, jac_triplets);
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
                normal_potential, normal_stiffness, dmin);
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
    const VectorMax2d u_aniso = collision.mu_aniso.cwiseProduct(u);

    return collision.weight * collision.normal_force_magnitude
        * mu_f0(u_aniso.norm(), collision.mu_s, collision.mu_k);
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
    const VectorMax2d u_aniso = collision.mu_aniso.cwiseProduct(u);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        collision.relative_velocity_matrix().transpose()
        * collision.tangent_basis;

    // Compute μ(‖u_aniso‖) f₁(‖u_aniso‖)/‖u_aniso‖
    const double mu_f1_over_norm_u =
        mu_f1_over_x(u_aniso.norm(), collision.mu_s, collision.mu_k);

    // μ(‖u_aniso‖) N(xᵗ) f₁(‖u_aniso‖)/‖u_aniso‖ T(xᵗ) u_aniso ∈ (n×2)(2×1) = (n×1)
    return T
        * ((collision.weight * collision.normal_force_magnitude
            * mu_f1_over_norm_u)
           * u_aniso);
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
    const VectorMax2d u_aniso = collision.mu_aniso.cwiseProduct(u);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        collision.relative_velocity_matrix().transpose()
        * collision.tangent_basis;

    // Compute ‖u_aniso‖
    const double norm_u = u_aniso.norm();

    // Compute μ(‖u_aniso‖) f₁(‖u_aniso‖)/‖u_aniso‖
    const double mu_f1_over_norm_u =
        mu_f1_over_x(norm_u, collision.mu_s, collision.mu_k);

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
            // I - u_aniso u_anisoᵀ/‖u_aniso‖² = ūūᵀ / ‖u_aniso‖² (where ū⋅u_aniso = 0)
            const Eigen::Vector2d u_perp(-u_aniso[1], u_aniso[0]);
            // Apply anisotropic scaling to T: T_aniso = T * diag(mu_aniso)
            MatrixMax<double, 12, 2> T_aniso = T;
            T_aniso.col(0) *= collision.mu_aniso[0];
            T_aniso.col(1) *= collision.mu_aniso[1];
            hess = // grouped to reduce number of operations
                (T_aniso * ((scale * mu_f1_over_norm_u / (norm_u * norm_u)) * u_perp))
                * (u_perp.transpose() * T_aniso.transpose());
        }
    } else if (norm_u == 0) {
        // ∇²D = μ N T [(f₂(‖u‖)‖u‖ − f₁(‖u‖))/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
        // lim_{‖u‖→0} ∇²D = μ N T [f₁(‖u‖)/‖u‖ I] Tᵀ
        // no PSD projection needed because μ N f₁(‖ū‖)/‖ū‖ ≥ 0
        if (project_hessian_to_psd != PSDProjectionMethod::NONE && scale <= 0) {
            hess.setZero(collision.ndof(), collision.ndof()); // -PSD = NSD ⟹ 0
        } else {
            // Apply anisotropic scaling to T: T_aniso = T * diag(mu_aniso)
            MatrixMax<double, 12, 2> T_aniso = T;
            T_aniso.col(0) *= collision.mu_aniso[0];
            T_aniso.col(1) *= collision.mu_aniso[1];
            hess = scale * mu_f1_over_norm_u * T_aniso * T_aniso.transpose();
        }
    } else {
        // ∇²D(v) = μ N T [f₂(‖u_aniso‖) u_aniso u_anisoᵀ + f₁(‖u_aniso‖)/‖u_aniso‖ I] Tᵀ
        //  ⟹ only need to project the inner 2x2 matrix to PSD
        const double f2 =
            mu_f2_x_minus_mu_f1_over_x3(norm_u, collision.mu_s, collision.mu_k);

        MatrixMax2d inner_hess = f2 * u_aniso * u_aniso.transpose();
        inner_hess.diagonal().array() += mu_f1_over_norm_u;
        inner_hess *= scale; // NOTE: negative scaling will be projected out
        inner_hess = project_to_psd(inner_hess, project_hessian_to_psd);

        // Apply anisotropic scaling to T: T_aniso = T * diag(mu_aniso)
        // This accounts for ∂u_aniso/∂u = diag(mu_aniso) in the chain rule
        MatrixMax<double, 12, 2> T_aniso = T;
        T_aniso.col(0) *= collision.mu_aniso[0];
        T_aniso.col(1) *= collision.mu_aniso[1];
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
    const double normal_stiffness,
    const double dmin,
    const bool no_mu) const
{
    // x is the rest position
    // u is the displacment at the begginging of the lagged solve
    // v is the current velocity
    //
    // τ = T(x + u)ᵀv is the tangential sliding velocity
    // τ_aniso = mu_aniso ⊙ τ is the anisotropically-scaled velocity
    // F(x, u, v) = -μ N(x + u) f₁(‖τ_aniso‖)/‖τ_aniso‖ T(x + u) τ_aniso
    //
    // Combined anisotropic friction model:
    //   1. mu_aniso velocity scaling: τ_aniso = diag(mu_aniso) · τ
    //   2. Direction-dependent coefficients (when mu_s_aniso.squaredNorm() > 0):
    //      - Direction computed from τ_aniso: τ_dir = τ_aniso / ‖τ_aniso‖
    //      - Effective mu from ellipse: μ_eff = ‖diag(μ_aniso) · τ_dir‖
    //   3. Isotropic path (when mu_s_aniso is zero): uses scalar mu_s/mu_k
    assert(rest_positions.size() == lagged_displacements.size());
    assert(rest_positions.size() == velocities.size());

    // const VectorMax12d x = dof(rest_positions, edges, faces);
    // const VectorMax12d u = dof(lagged_displacements, edges, faces);
    // const VectorMax12d v = dof(velocities, edges, faces);

    // x:
    const VectorMax12d lagged_positions = rest_positions + lagged_displacements;

    // Compute N(x + u)
    const double N = normal_potential.force_magnitude(
        collision.compute_distance(lagged_positions), dmin, normal_stiffness);

    // Compute P
    const MatrixMax<double, 3, 2> P =
        collision.compute_tangent_basis(lagged_positions);

    // compute β
    const VectorMax2d beta = collision.compute_closest_point(lagged_positions);

    // Compute Γ
    const MatrixMax<double, 3, 12> Gamma =
        collision.relative_velocity_matrix(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;

    // Compute τ = PᵀΓv
    const VectorMax2d tau = T.transpose() * velocities;

    // Always apply mu_aniso velocity scaling first: tau_aniso = mu_aniso ⊙ tau
    const VectorMax2d tau_aniso = collision.mu_aniso.cwiseProduct(tau);

    // Compute effective mu (handles both anisotropic and isotropic cases)
    const auto [mu_s, mu_k] = compute_anisotropic_mu_eff_from_tau_aniso(
        tau_aniso, collision.mu_s_aniso, collision.mu_k_aniso,
        collision.mu_s, collision.mu_k, no_mu);

    // Compute μ(‖τ_aniso‖) f₁(‖τ_aniso‖)/‖τ_aniso‖
    const double tau_aniso_norm = tau_aniso.norm();
    const double mu_f1_over_norm_tau = mu_f1_over_x(tau_aniso_norm, mu_s, mu_k);

    // F = -μ N f₁(‖τ_aniso‖)/‖τ_aniso‖ T τ_aniso
    // NOTE: no_mu -> leave mu out of this function (i.e., assuming mu = 1)
    // NOTE: Force always uses tau_aniso (with mu_aniso scaling applied).
    //       When anisotropic friction is enabled, mu_eff is computed from
    //       the direction of tau_aniso.
    return -collision.weight * N * mu_f1_over_norm_tau * T * tau_aniso;
}

MatrixMax12d TangentialPotential::force_jacobian(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMax12d> rest_positions,       // = x
    Eigen::ConstRef<VectorMax12d> lagged_displacements, // = u
    Eigen::ConstRef<VectorMax12d> velocities,           // = v
    const NormalPotential& normal_potential,
    const double normal_stiffness,
    const DiffWRT wrt,
    const double dmin) const
{
    // x is the rest position
    // u is the displacment at the begginging of the lagged solve
    // v is the current velocity
    //
    // τ = T(x + u)ᵀv is the tangential sliding velocity
    // τ_aniso = mu_aniso ⊙ τ is the anisotropically-scaled velocity
    // F(x, u, v) = -μ N(x + u) f₁(‖τ_aniso‖)/‖τ_aniso‖ T(x + u) τ_aniso
    //
    // Compute ∇F with combined anisotropic friction model:
    //   - mu_aniso velocity scaling is ALWAYS applied via jac_tau_aniso
    //   - When mu_s_aniso.squaredNorm() > 0: effective mu depends on
    //     tau_aniso direction, so d(mu_eff)/d(tau_aniso) is included
    //   - When mu_s_aniso is zero (isotropic): uses scalar mu_s/mu_k
    assert(rest_positions.size() == lagged_displacements.size());
    assert(lagged_displacements.size() == velocities.size());
    const int n = rest_positions.size();
    const int dim = n / collision.num_vertices();
    assert(n % collision.num_vertices() == 0);

    // const VectorMax12d x = dof(rest_positions, edges, faces);
    // const VectorMax12d u = dof(lagged_displacements, edges, faces);
    // const VectorMax12d v = dof(velocities, edges, faces);

    // x + u:
    const VectorMax12d lagged_positions = rest_positions + lagged_displacements;
    const bool need_jac_N_or_T = wrt != DiffWRT::VELOCITIES;

    // Compute N
    const double N = normal_potential.force_magnitude(
        collision.compute_distance(lagged_positions), dmin, normal_stiffness);

    // Compute ∇N
    VectorMax12d grad_N;
    if (need_jac_N_or_T) {
        // ∇ₓN = ∇ᵤN
        grad_N = normal_potential.force_magnitude_gradient(
            collision.compute_distance(lagged_positions),
            collision.compute_distance_gradient(lagged_positions), dmin,
            normal_stiffness);
        assert(grad_N.array().isFinite().all());
    }

    // Compute P
    const MatrixMax<double, 3, 2> P =
        collision.compute_tangent_basis(lagged_positions);

    // Compute β
    const VectorMax2d beta = collision.compute_closest_point(lagged_positions);

    // Compute Γ
    const MatrixMax<double, 3, 12> Gamma =
        collision.relative_velocity_matrix(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T = Gamma.transpose() * P;

    // Compute ∇T
    MatrixMax<double, 144, 2> jac_T;
    if (need_jac_N_or_T) {
        jac_T.resize(n * n, dim - 1);
        // ∇T = ∇(ΓᵀP) = ∇ΓᵀP + Γᵀ∇P
        const MatrixMax<double, 36, 2> jac_P =
            collision.compute_tangent_basis_jacobian(lagged_positions);
        for (int i = 0; i < n; i++) {
            // ∂T/∂xᵢ += Γᵀ ∂P/∂xᵢ
            jac_T.middleRows(i * n, n) =
                Gamma.transpose() * jac_P.middleRows(i * dim, dim);
        }

        // Vertex-vertex does not have a closest point
        if (beta.size() != 0) {
            // ∇Γ(β) = ∇ᵦΓ∇β ∈ ℝ^{d×n×n} ≡ ℝ^{nd×n}
            const MatrixMax<double, 2, 12> jac_beta =
                collision.compute_closest_point_jacobian(lagged_positions);
            const MatrixMax<double, 6, 12> jac_Gamma_wrt_beta =
                collision.relative_velocity_matrix_jacobian(beta);

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
    const VectorMax2d tau = P.transpose() * Gamma * velocities;

    // Always apply mu_aniso velocity scaling first: tau_aniso = mu_aniso ⊙ tau
    const VectorMax2d tau_aniso = collision.mu_aniso.cwiseProduct(tau);

    // Compute ∇τ = ∇T(x + u)ᵀv + T(x + u)ᵀ∇v
    MatrixMax<double, 2, 12> jac_tau;
    if (need_jac_N_or_T) {
        jac_tau.resize(dim - 1, n);
        // Compute ∇T(x + u)ᵀv
        for (int i = 0; i < n; i++) {
            jac_tau.col(i) =
                jac_T.middleRows(i * n, n).transpose() * velocities;
        }
    } else {
        jac_tau = T.transpose(); // Tᵀ ∇ᵥv = Tᵀ
    }

    // Compute ∇tau_aniso = diag(mu_aniso) * ∇tau (chain rule for mu_aniso
    // scaling)
    MatrixMax<double, 2, 12> jac_tau_aniso = jac_tau;
    jac_tau_aniso.row(0) *= collision.mu_aniso[0];
    jac_tau_aniso.row(1) *= collision.mu_aniso[1];

    // Check if direction-dependent anisotropic friction is enabled
    const bool is_anisotropic = collision.mu_s_aniso.squaredNorm() > 0;

    // Compute effective mu (handles both anisotropic and isotropic cases)
    const auto [mu_s, mu_k] = compute_anisotropic_mu_eff_from_tau_aniso(
        tau_aniso, collision.mu_s_aniso, collision.mu_k_aniso,
        collision.mu_s, collision.mu_k, false);

    // Compute μ f₁(‖τ_aniso‖)/‖τ_aniso‖
    const double tau_aniso_norm = tau_aniso.norm();
    const double mu_f1_over_norm_tau = mu_f1_over_x(tau_aniso_norm, mu_s, mu_k);

    // Compute ∇(μ f₁(‖τ_aniso‖)/‖τ_aniso‖)
    VectorMax12d grad_mu_f1_over_norm_tau;
    if (tau_aniso_norm == 0) {
        // lim_{x→0} f₂(x)x² = 0
        grad_mu_f1_over_norm_tau.setZero(n);
    } else {
        if (is_anisotropic) {
            // For direction-dependent anisotropic friction, we need to include
            // d(mu_eff)/d(tau_aniso) term. The combined model computes
            // direction from tau_aniso.

            // Main term: (f₂(‖τ_aniso‖)‖τ_aniso‖ - f₁(‖τ_aniso‖)) /
            // ‖τ_aniso‖³ τ_anisoᵀ ∇τ_aniso. This treats mu_eff as constant
            // (evaluated at current tau_aniso direction).
            double f2 =
                mu_f2_x_minus_mu_f1_over_x3(tau_aniso_norm, mu_s, mu_k);
            assert(std::isfinite(f2));
            grad_mu_f1_over_norm_tau =
                f2 * tau_aniso.transpose() * jac_tau_aniso;

            // Additional term: d(mu_eff)/d(tau_aniso) contribution.
            // For the combined model, the derivative is computed with respect
            // to tau_aniso (the scaled velocity), not raw tau.
            const auto [dmu_s_eff_dtau, dmu_k_eff_dtau] =
                compute_anisotropic_mu_eff_derivatives(
                    tau_aniso, collision.mu_s_aniso, collision.mu_k_aniso,
                    mu_s, mu_k);

            // Approximate the contribution: ∂(μ f₁/‖τ‖)/∂μ_eff * dμ_eff/dτ *
            // ∇τ. We use the average of static and kinetic derivatives as an
            // approximation.
            Eigen::Vector2d dmu_eff_dtau_avg =
                (dmu_s_eff_dtau + dmu_k_eff_dtau) * 0.5;

            // The derivative of mu_f1_over_x with respect to mu_eff is
            // approximately smooth_mu_f1_over_x evaluated at the current
            // point. We multiply by the change in mu_eff per unit change in
            // tau_aniso direction.
            VectorMax12d dmu_eff_contribution =
                dmu_eff_dtau_avg.transpose() * jac_tau_aniso;
            // Scale by the sensitivity of mu_f1_over_x to changes in mu
            grad_mu_f1_over_norm_tau +=
                0.1 * mu_f1_over_norm_tau * dmu_eff_contribution;
        } else {
            // Isotropic: ∇ (f₁(‖tau_aniso‖)/‖tau_aniso‖) = (f₂‖τ‖ - f₁) / ‖τ‖³
            // τ_anisoᵀ ∇tau_aniso
            double f2 =
                mu_f2_x_minus_mu_f1_over_x3(tau_aniso_norm, mu_s, mu_k);
            assert(std::isfinite(f2));
            grad_mu_f1_over_norm_tau =
                f2 * tau_aniso.transpose() * jac_tau_aniso;
        }
    }

    // Premultiplied values
    const VectorMax12d T_times_tau = T * tau_aniso;

    // ------------------------------------------------------------------------
    // Compute J = ∇F = ∇(-μ N f₁(‖τ_aniso‖)/‖τ_aniso‖ T τ_aniso)
    MatrixMax12d J = MatrixMax12d::Zero(n, n);

    // = -μ f₁(‖τ_aniso‖)/‖τ_aniso‖ (T τ_aniso) [∇N]ᵀ
    if (need_jac_N_or_T) {
        J = mu_f1_over_norm_tau * T_times_tau * grad_N.transpose();
    }

    // + -N T τ_aniso [∇(μ f₁(‖τ_aniso‖)/‖τ_aniso‖)]
    J += N * T_times_tau * grad_mu_f1_over_norm_tau.transpose();

    // + -μ N f₁(‖τ_aniso‖)/‖τ_aniso‖ [∇T] τ_aniso
    if (need_jac_N_or_T) {
        const VectorMax2d scaled_tau = N * mu_f1_over_norm_tau * tau_aniso;
        for (int i = 0; i < n; i++) {
            // ∂J/∂xᵢ = ∂T/∂xᵢ * τ_aniso
            J.col(i) += jac_T.middleRows(i * n, n) * scaled_tau;
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
                // smooth_contact_force_jacobian(
                //     collision, collision.dof(rest_positions, edges, faces),
                //     collision.dof(lagged_displacements, edges, faces),
                //     collision.dof(velocities, edges, faces), //
                //     wrt);

                const auto vis =
                    collision.vertex_ids(mesh.edges(), mesh.faces());

                local_hessian_to_global_triplets(
                    local_force_jacobian, vis, dim, jac_triplets);

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
                    cc_vert_ids, dim, jac_triplets);
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

    // This term is zero!
    // if wrt == X then compute ∇ₓ w(x)
    // if (wrt == DiffWRT::REST_POSITIONS) {
    //     for (int i = 0; i < collisions.size(); i++) {
    //         const FrictionCollision& collision = collisions[i];
    //         assert(collision.weight_gradient.size() ==
    //         rest_positions.size()); if (collision.weight_gradient.size() !=
    //         rest_positions.size()) {
    //             throw std::runtime_error(
    //                 "Shape derivative is not computed for friction
    //                 collision!");
    //         }

    //         VectorMaxNd local_force
    //         =
    //             smooth_contact_force(
    //                 collision, collision.dof(rest_positions, edges, faces),
    //                 collision.dof(lagged_displacements, edges, faces),
    //                 collision.dof(velocities, edges, faces), //
    //                 dmin);
    //         assert(collision.weight != 0);
    //         local_force /= collision.weight;

    //         Eigen::SparseVector<double> force(rest_positions.size());
    //         force.reserve(local_force.size());
    //         local_gradient_to_global_gradient(
    //             local_force, collision.vertex_ids(edges, faces), dim, force);

    //         jacobian += force * collision.weight_gradient.transpose();
    //     }
    // }

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
    //     lagged_positions, barrier_potential, barrier_stiffness, dmin);
    const double N = collision.normal_force_magnitude;

    // Compute P
    const MatrixMax<double, 3, 2> P =
        collision.compute_tangent_basis(lagged_positions);

    // compute β
    const VectorMax2d beta = collision.compute_closest_point(lagged_positions);

    // Compute Γ
    const MatrixMax<double, 3, STENCIL_NDOF> Gamma =
        collision.relative_velocity_matrix(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, STENCIL_NDOF, 2> T = Gamma.transpose() * P;

    // Compute τ = PᵀΓv
    const VectorMax2d tau = T.transpose() * velocities;

    // Apply anisotropic scaling: tau_aniso = mu_aniso ⊙ tau
    const VectorMax2d tau_aniso = collision.mu_aniso.cwiseProduct(tau);

    // Compute f₁(‖tau_aniso‖)/‖tau_aniso‖
    const double mu_s = no_mu ? 1.0 : collision.mu_s;
    const double mu_k = no_mu ? 1.0 : collision.mu_k;
    const double mu_f1_over_norm_tau = mu_f1_over_x(tau_aniso.norm(), mu_s, mu_k);

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
    const int n = lagged_positions.size();
    const int dim = n / collision.num_vertices();
    assert(n % collision.num_vertices() == 0);

    const bool need_jac_N_or_T = wrt != DiffWRT::VELOCITIES;

    // Compute P
    const MatrixMax<double, 3, 2> P =
        collision.compute_tangent_basis(lagged_positions);

    // Compute β
    const VectorMax2d beta = collision.compute_closest_point(lagged_positions);

    // Compute Γ
    const MatrixMax<double, 3, STENCIL_NDOF> Gamma =
        collision.relative_velocity_matrix(beta);

    // Compute T = ΓᵀP
    const MatrixMax<double, STENCIL_NDOF, 2> T = Gamma.transpose() * P;

    // Compute ∇T
    MatrixMax<double, STENCIL_NDOF * STENCIL_NDOF, 2> jac_T;
    if (need_jac_N_or_T) {
        jac_T.resize(n * n, dim - 1);
        // ∇T = ∇(ΓᵀP) = ∇ΓᵀP + Γᵀ∇P
        const MatrixMax<double, 36, 2> jac_P =
            collision.compute_tangent_basis_jacobian(lagged_positions);
        for (int i = 0; i < n; i++) {
            // ∂T/∂xᵢ += Γᵀ ∂P/∂xᵢ
            jac_T.middleRows(i * n, n) =
                Gamma.transpose() * jac_P.middleRows(i * dim, dim);
        }

        // Vertex-vertex does not have a closest point
        if (beta.size()) {
            // ∇Γ(β) = ∇ᵦΓ∇β ∈ ℝ^{d×n×n} ≡ ℝ^{nd×n}
            const MatrixMax<double, 2, STENCIL_NDOF> jac_beta =
                collision.compute_closest_point_jacobian(lagged_positions);
            const MatrixMax<double, 6, STENCIL_NDOF> jac_Gamma_wrt_beta =
                collision.relative_velocity_matrix_jacobian(beta);

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
    const VectorMax2d tau = P.transpose() * Gamma * velocities;

    // Apply anisotropic scaling: tau_aniso = mu_aniso ⊙ tau
    const VectorMax2d tau_aniso = collision.mu_aniso.cwiseProduct(tau);

    // Compute ∇τ = ∇T(x + u)ᵀv + T(x + u)ᵀ∇v
    MatrixMax<double, 2, STENCIL_NDOF> jac_tau;
    if (need_jac_N_or_T) {
        jac_tau.resize(dim - 1, n);
        // Compute ∇T(x + u)ᵀv
        for (int i = 0; i < n; i++) {
            jac_tau.col(i) =
                jac_T.middleRows(i * n, n).transpose() * velocities;
        }
    } else {
        jac_tau = T.transpose(); // Tᵀ ∇ᵥv = Tᵀ
    }

    // Compute ∇tau_aniso = diag(mu_aniso) * ∇tau
    MatrixMax<double, 2, STENCIL_NDOF> jac_tau_aniso = jac_tau;
    jac_tau_aniso.row(0) *= collision.mu_aniso[0];
    jac_tau_aniso.row(1) *= collision.mu_aniso[1];

    // Compute f₁(‖tau_aniso‖)/‖tau_aniso‖
    const double tau_norm = tau_aniso.norm();
    const double mu_s = no_mu ? 1.0 : collision.mu_s;
    const double mu_k = no_mu ? 1.0 : collision.mu_k;
    const double f1_over_norm_tau = mu_f1_over_x(tau_norm, mu_s, mu_k);

    // Compute ∇(f₁(‖tau_aniso‖)/‖tau_aniso‖)
    VectorMaxNd grad_f1_over_norm_tau;
    if (tau_norm == 0) {
        // lim_{x→0} f₂(x)x² = 0
        grad_f1_over_norm_tau.setZero(n);
    } else {
        // ∇ (f₁(‖tau_aniso‖)/‖tau_aniso‖) = (f₁'(‖tau_aniso‖)‖tau_aniso‖ - f₁(‖tau_aniso‖)) / ‖tau_aniso‖³ tau_anisoᵀ ∇tau_aniso
        double f2 = mu_f2_x_minus_mu_f1_over_x3(tau_norm, mu_s, mu_k);
        assert(std::isfinite(f2));
        grad_f1_over_norm_tau = f2 * tau_aniso.transpose() * jac_tau_aniso;
    }

    // Premultiplied values
    const VectorMaxNd T_times_tau = T * tau_aniso;

    // ------------------------------------------------------------------------
    // Compute J = ∇F = ∇(-μ N f₁(‖τ‖)/‖τ‖ T τ)
    MatrixMaxNd J = MatrixMaxNd::Zero(n, n);

    // + -μ N T τ [∇(f₁(‖τ‖)/‖τ‖)]
    J += T_times_tau * grad_f1_over_norm_tau.transpose();

    // + -μ N f₁(‖tau_aniso‖)/‖tau_aniso‖ [∇T] tau_aniso
    if (need_jac_N_or_T) {
        const VectorMax2d scaled_tau = f1_over_norm_tau * tau_aniso;
        for (int i = 0; i < n; i++) {
            // ∂J/∂xᵢ = ∂T/∂xᵢ * tau_aniso
            J.col(i) += jac_T.middleRows(i * n, n) * scaled_tau;
        }
    }

    // + -μ N f₁(‖tau_aniso‖)/‖tau_aniso‖ T [∇tau_aniso]
    J += f1_over_norm_tau * T * jac_tau_aniso;

    // NOTE: ∇ₓw(x) is not local to the collision pair (i.e., it involves more
    // than the 4 collisioning vertices), so we do not have enough information
    // here to compute the gradient. Instead this should be handled outside of
    // the function. For a simple multiplicitive model (∑ᵢ wᵢ Fᵢ) this can be
    // done easily.
    J *= -collision.weight;

    return J;
}

} // namespace ipc