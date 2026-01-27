#include <common.hpp>

#include <ipc/friction/smooth_mu.hpp>

using namespace ipc;

void define_smooth_mu(py::module_& m)
{
    m.def(
        "smooth_mu", &smooth_mu,
        R"ipc_Qu8mg5v7(
        Smooth coefficient from static to kinetic friction.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static friction.
            mu_k: Coefficient of kinetic friction.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the μ at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_v"_a);

    m.def(
        "smooth_mu_derivative", &smooth_mu_derivative,
        R"ipc_Qu8mg5v7(
        Compute the derivative of the smooth coefficient from static to kinetic friction.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static friction.
            mu_k: Coefficient of kinetic friction.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the derivative at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_v"_a);

    m.def(
        "smooth_mu_f0", &smooth_mu_f0,
        R"ipc_Qu8mg5v7(
        Compute the value of the ∫ μ(y) f₁(y) dy, where f₁ is the first derivative of the smooth friction mollifier.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static friction.
            mu_k: Coefficient of kinetic friction.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the integral at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_v"_a);

    m.def(
        "smooth_mu_f1", &smooth_mu_f1,
        R"ipc_Qu8mg5v7(
        Compute the value of the μ(y) f₁(y), where f₁ is the first derivative of the smooth friction mollifier.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static friction.
            mu_k: Coefficient of kinetic friction.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the product at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_v"_a);

    m.def(
        "smooth_mu_f2", &smooth_mu_f2,
        R"ipc_Qu8mg5v7(
        Compute the value of d/dy (μ(y) f₁(y)), where f₁ is the first derivative of the smooth friction mollifier.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static friction.
            mu_k: Coefficient of kinetic friction.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the derivative at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_v"_a);

    m.def(
        "smooth_mu_f1_over_x", &smooth_mu_f1_over_x,
        R"ipc_Qu8mg5v7(
        Compute the value of the μ(y) f₁(y) / y, where f₁ is the first derivative of the smooth friction mollifier.

        Note:
            The `x` in the function name refers to the parameter `y`.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static friction.
            mu_k: Coefficient of kinetic friction.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the product at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_v"_a);

    m.def(
        "smooth_mu_f2_x_minus_mu_f1_over_x3",
        &smooth_mu_f2_x_minus_mu_f1_over_x3,
        R"ipc_Qu8mg5v7(
        Compute the value of the [(d/dy μ(y) f₁(y)) ⋅ y - μ(y) f₁(y)] / y³, where f₁ and f₂ are the first and second derivatives of the smooth friction mollifier.

        Note:
            The `x` in the function name refers to the parameter `y`.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static friction.
            mu_k: Coefficient of kinetic friction.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the expression at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_v"_a);

    m.def(
        "anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq",
        &anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq,
        R"ipc_Qu8mg5v7(
        Compute effective friction coefficients for elliptical anisotropy (L2
        projection).

        For anisotropic friction, the friction coefficient depends on the
        direction of tangential velocity. This function computes the effective
        friction coefficients along a given direction using the elliptical L2
        projection model: μ_eff = sqrt((μ₀ t₀)² + (μ₁ t₁)²), where t = τ / ||τ||
        is the unit direction vector. The function name describes the mathematical
        expression: sqrt(mu0*t0² + mu1*t1²).

        Parameters:
            tau_dir: Unit 2D direction of tangential velocity (tau / ||tau||).
                     Must be a unit vector.
            mu_s_aniso: Static friction ellipse axes (2D vector). Each component
                        represents the friction coefficient along the
                        corresponding tangent basis direction.
            mu_k_aniso: Kinetic friction ellipse axes (2D vector). Each
                        component represents the friction coefficient along the
                        corresponding tangent basis direction.

        Returns:
            A tuple containing (mu_s_eff, mu_k_eff), the effective static and
            kinetic friction coefficients along the direction tau_dir.

        Note:
            If mu_s_aniso and mu_k_aniso are zero vectors, the function returns
            (0, 0), which triggers backward-compatible isotropic behavior.
        )ipc_Qu8mg5v7",
        "tau_dir"_a, "mu_s_aniso"_a, "mu_k_aniso"_a);

    m.def(
        "anisotropic_mu_eff_dtau", &anisotropic_mu_eff_dtau,
        R"ipc_Qu8mg5v7(
        Compute the derivative of effective friction coefficient with respect
        to tangential velocity for elliptical anisotropy.

        This function computes ∂μ_eff/∂τ for the elliptical anisotropy model.
        The derivative is needed for computing the Jacobian of friction forces
        when anisotropic friction is enabled.

        Parameters:
            tau: Tangential velocity (2D vector). The velocity vector in the
                 tangent plane.
            mu_aniso: Ellipse axes (2D vector). The anisotropic friction
                      coefficients along each tangent direction.
            mu_eff: Effective friction coefficient computed from
                    anisotropic_mu_eff_sqrt_mu0_t0_sq_plus_mu1_t1_sq(). This is
                    passed to avoid recomputation.

        Returns:
            The derivative ∂μ_eff/∂τ as a 2D vector.

        Note:
            Returns zero vector if ||tau|| ≈ 0 or mu_eff ≈ 0 to handle edge
            cases gracefully.
        )ipc_Qu8mg5v7",
        "tau"_a, "mu_aniso"_a, "mu_eff"_a);
}
