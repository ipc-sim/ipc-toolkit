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
        "anisotropic_mu_eff_f", &anisotropic_mu_eff_f,
        R"ipc_Qu8mg5v7(
        Effective static and kinetic friction along a unit direction for the
        elliptical (matchstick) model: μ_eff = sqrt((μ₀ t₀)² + (μ₁ t₁)²).
        Matchstick model: Erleben et al., CGF 2019, DOI 10.1111/cgf.13885.

        Parameters:
            tau_dir: Unit 2D direction (tau / ||tau||).
            mu_s_aniso: Static friction ellipse axes (2D).
            mu_k_aniso: Kinetic friction ellipse axes (2D).

        Returns:
            (mu_s_eff, mu_k_eff) along tau_dir. (0, 0) if inputs are zero
            (isotropic fallback).
        )ipc_Qu8mg5v7",
        "tau_dir"_a, "mu_s_aniso"_a, "mu_k_aniso"_a);

    m.def(
        "anisotropic_mu_eff_f_dtau", &anisotropic_mu_eff_f_dtau,
        R"ipc_Qu8mg5v7(
        ∂μ_eff/∂τ for the elliptical model (friction force Jacobians).
        Matchstick model: Erleben et al., CGF 2019, DOI 10.1111/cgf.13885.

        Parameters:
            tau: Tangential velocity (2D) in the tangent plane.
            mu_aniso: Ellipse axes (2D).
            mu_eff: Effective μ from anisotropic_mu_eff_f (avoids recomputation).

        Returns:
            ∂μ_eff/∂τ as 2D vector. Zero if ||tau|| ≈ 0 or mu_eff ≈ 0.
        )ipc_Qu8mg5v7",
        "tau"_a, "mu_aniso"_a, "mu_eff"_a);
}
