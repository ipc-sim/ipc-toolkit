#include "friction_potential.hpp"
#include <ipc/tangent/tangent_basis.hpp>  // Add this include for TangentBasis

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

                const std::array<long, 4> vis =
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

                const std::array<long, 4> vis =
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
    // μ N(xᵗ) f₀(‖u‖) (where u = T(xᵗ)ᵀv)

    // Compute u = PᵀΓv
    const VectorMax2d u = collision.tangent_basis.transpose()
        * collision.relative_velocity(velocities);

    return collision.weight * collision.mu * collision.normal_force_magnitude *
        ((collision.s_mu > 0 && collision.k_mu > 0) ? 
            f0_mus(u.norm(), collision.s_mu, collision.k_mu) : 
            f0(u.norm()));
}

VectorMax12d TangentialPotential::gradient(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMax12d> velocities) const
{
    // ∇ₓ μ N(xᵗ) f₀(‖u‖) (where u = T(xᵗ)ᵀv)
    //  = μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u

    // Compute u = PᵀΓv
    const VectorMax2d u = collision.tangent_basis.transpose()
        * collision.relative_velocity(velocities);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        collision.relative_velocity_matrix().transpose()
        * collision.tangent_basis;

    // Compute f₁(‖u‖)/‖u‖
    const double f1_over_norm_u = (collision.s_mu > 0 && collision.k_mu > 0) ?
                                f1_over_x_mus(u.norm(), collision.s_mu, collision.k_mu) :
                                f1_over_x(u.norm());

    // μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u ∈ ℝⁿ
    // (n×2)(2×1) = (n×1)
    double mu = (is_dynamic(u.norm()) && collision.k_mu != -1)
                    ? collision.k_mu
                    : (collision.s_mu != -1 ? collision.s_mu : collision.mu);

    return T
        * ((collision.weight * mu * collision.normal_force_magnitude
            * f1_over_norm_u)
           * u);
}

Eigen::VectorXd TangentialPotential::gradient(
    const Eigen::MatrixXd& V,
    const Eigen::VectorXd& u0,
    const Eigen::Ref<const Eigen::MatrixXd>& basis,  // Updated parameter type
    double dt,
    double mu,
    double eps_v,
    double static_mu, // Added parameter for static friction coefficient
    double kinetic_mu) const // Added parameter for kinetic friction coefficient
{
    // If static and kinetic mu are not specified (or are the same),
    // use the original implementation
    if (static_mu <= 0 || std::abs(static_mu - kinetic_mu) < 1e-10) {
        // Original implementation for single mu
        const Eigen::VectorXd u0_local = basis.transpose() * u0;
        const Eigen::VectorXd u = u0_local;
        const double u_norm = u.norm();

        if (u_norm < 1e-10) {
            return Eigen::VectorXd::Zero(u0.size());
        }

        const double force_scale = mu * f1_over_x(u_norm) / u_norm;
        return basis * (force_scale * u);
    }
    
    // Use the smooth transition between static and dynamic friction
    const Eigen::VectorXd u0_local = basis.transpose() * u0;
    const Eigen::VectorXd u = u0_local;
    const double u_norm = u.norm();

    if (u_norm < 1e-10) {
        return Eigen::VectorXd::Zero(u0.size());
    }

    // Simple transition from static to kinetic friction
    double mu_effective = kinetic_mu;
    if (u_norm < eps_v) {
        // Blend between static and kinetic friction based on velocity
        double t = u_norm / eps_v; // 0 to 1
        mu_effective = static_mu + (kinetic_mu - static_mu) * t;
    }
    
    const double force_scale = mu_effective * f1_over_x(u_norm) / u_norm;
    return basis * (force_scale * u);
}

MatrixMax12d TangentialPotential::hessian(
    const TangentialCollision& collision,
    Eigen::ConstRef<VectorMax12d> velocities,
    const PSDProjectionMethod project_hessian_to_psd) const
{
    // ∇ₓ μ N(xᵗ) f₁(‖u‖)/‖u‖ T(xᵗ) u (where u = T(xᵗ)ᵀ v)
    //  = μ N T [(f₂(‖u‖)‖u‖ − f₁(‖u‖))/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
    //  = μ N T [f₂(‖u‖) uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ

    // Compute u = PᵀΓv
    const VectorMax2d u = collision.tangent_basis.transpose()
        * collision.relative_velocity(velocities);

    // Compute T = ΓᵀP
    const MatrixMax<double, 12, 2> T =
        collision.relative_velocity_matrix().transpose()
        * collision.tangent_basis;

    // Compute ‖u‖
    const double norm_u = u.norm();

    // Compute f₁(‖u‖)/‖u‖
    const double f1_over_norm_u = (collision.s_mu > 0 && collision.k_mu > 0) ?
                                f1_over_x_mus(norm_u, collision.s_mu, collision.k_mu) :
                                f1_over_x(norm_u);

    // Compute μ N(xᵗ)
    double mu = (is_dynamic(u.norm()) && collision.k_mu != -1)
                    ? collision.k_mu
                    : (collision.s_mu != -1 ? collision.s_mu : collision.mu);

    const double scale =
        collision.weight * mu * collision.normal_force_magnitude;

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
            // I - uuᵀ/‖u‖² = ūūᵀ / ‖u‖² (where ū⋅u = 0)
            const Eigen::Vector2d u_perp(-u[1], u[0]);
            hess = // grouped to reduce number of operations
                (T * ((scale * f1_over_norm_u / (norm_u * norm_u)) * u_perp))
                * (u_perp.transpose() * T.transpose());
        }
    } else if (norm_u == 0) {
        // ∇²D = μ N T [(f₂(‖u‖)‖u‖ − f₁(‖u‖))/‖u‖³ uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
        // lim_{‖u‖→0} ∇²D = μ N T [f₁(‖u‖)/‖u‖ I] Tᵀ
        // no PSD projection needed because μ N f₁(‖ū‖)/‖ū‖ ≥ 0
        if (project_hessian_to_psd != PSDProjectionMethod::NONE && scale <= 0) {
            hess.setZero(collision.ndof(), collision.ndof()); // -PSD = NSD ⟹ 0
        } else {
            hess = scale * f1_over_norm_u * T * T.transpose();
        }
    } else {
        // ∇²D(v) = μ N T [f₂(‖u‖) uuᵀ + f₁(‖u‖)/‖u‖ I] Tᵀ
        //  ⟹ only need to project the inner 2x2 matrix to PSD
        double f2 = (collision.s_mu > 0 && collision.k_mu > 0) ?
            f2_x_minus_f1_over_x3_mus(norm_u, collision.s_mu, collision.k_mu) :
            f2_x_minus_f1_over_x3(norm_u);

        MatrixMax2d inner_hess = f2 * u * u.transpose();
        inner_hess.diagonal().array() += f1_over_norm_u;
        inner_hess *= scale; // NOTE: negative scaling will be projected out
        inner_hess = project_to_psd(inner_hess, project_hessian_to_psd);

        hess = T * inner_hess * T.transpose();
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
    // F(x, u, v) = -μ N(x + u) f₁(‖τ‖)/‖τ‖ T(x + u) τ
    assert(rest_positions.size() == lagged_displacements.size());
    assert(rest_positions.size() == velocities.size());

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

    // Compute f₁(‖τ‖)/‖τ‖
    const double tau_norm = tau.norm();
    
    // Determine which mu to use
    double mu = collision.mu;
    
    // Function to look up material-specific friction coefficient
    static thread_local std::shared_ptr<Eigen::MatrixXd> material_friction_table = nullptr;
    
    if (!no_mu) {
        // For material-specific friction, look up the appropriate mu from material IDs
        if (collision.has_material_ids() && material_friction_table) {
            mu = get_friction_coefficient(collision, *material_friction_table);
        } else if (is_dynamic(tau_norm) && collision.k_mu > 0) {
            // Use kinetic friction if we're in dynamic regime
            mu = collision.k_mu;
        } else if (collision.s_mu > 0) {
            // Use static friction if we're in static regime
            mu = collision.s_mu;
        }
    } else {
        mu = 1.0; // no_mu case
    }
    
    // Calculate force mollifier factor
    double f1_over_norm_tau;
    if (collision.s_mu > 0 && collision.k_mu > 0) {
        f1_over_norm_tau = f1_over_x_mus(tau_norm, collision.s_mu, collision.k_mu);
    } else {
        f1_over_norm_tau = f1_over_x(tau_norm);
    }
    
    // F = -μ N f₁(‖τ‖)/‖τ‖ T τ
    return -collision.weight * mu * N * f1_over_norm_tau * T * tau;
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
    // F(x, u, v) = -μ N(x + u) f₁(‖τ‖)/‖τ‖ T(x + u) τ
    //
    // Compute ∇F
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
        if (beta.size()) {
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

    // Compute f₁(‖τ‖)/‖τ‖
    const double tau_norm = tau.norm();
    double f1_over_norm_tau = (collision.s_mu > 0 && collision.k_mu > 0) ?
                            f1_over_x_mus(tau_norm, collision.s_mu, collision.k_mu) :
                            f1_over_x(tau_norm);

    // Compute ∇(f₁(‖τ‖)/‖τ‖)
    VectorMax12d grad_f1_over_norm_tau;
    if (tau_norm == 0) {
        // lim_{x→0} f₂(x)x² = 0
        grad_f1_over_norm_tau.setZero(n);
    } else {
        // ∇ (f₁(‖τ‖)/‖τ‖) = (f₂(‖τ‖)‖τ‖ - f₁(‖τ‖)) / ‖τ‖³ τᵀ ∇τ
        double f2 = (collision.s_mu > 0 && collision.k_mu > 0) ?
            f2_x_minus_f1_over_x3_mus(tau_norm, collision.s_mu, collision.k_mu) :
            f2_x_minus_f1_over_x3(tau_norm);
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

    // NOTE: ∇ₓw(x) is not local to the collision pair (i.e., it involves more
    // than the 4 collisioning vertices), so we do not have enough information
    // here to compute the gradient. Instead this should be handled outside of
    // the function. For a simple multiplicitive model (∑ᵢ wᵢ Fᵢ) this can be
    // done easily.
    J *= -collision.weight * collision.mu;

    return J;
}

double TangentialPotential::get_friction_coefficient(
    const TangentialCollision& collision,
    const Eigen::MatrixXd& material_friction_table) const
{
    // If material IDs are not set, use the collision's mu value
    if (collision.material_id1 == NO_MATERIAL_ID || 
        collision.material_id2 == NO_MATERIAL_ID || 
        material_friction_table.size() == 0) {
        return collision.mu;
    }

    // Check that material IDs are valid indices for the friction table
    if (collision.material_id1 >= 0 && 
        collision.material_id1 < material_friction_table.rows() &&
        collision.material_id2 >= 0 && 
        collision.material_id2 < material_friction_table.cols()) {
        return material_friction_table(collision.material_id1, collision.material_id2);
    }
    
    // Default to collision's mu if material IDs are out of range
    return collision.mu;
}

} // namespace ipc