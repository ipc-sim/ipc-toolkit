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
}
