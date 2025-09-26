#include <common.hpp>

#include <ipc/adhesion/adhesion.hpp>

using namespace ipc;

void define_adhesion(py::module_& m)
{
    m.def(
        "normal_adhesion_potential", &normal_adhesion_potential,
        R"ipc_Qu8mg5v7(
        The normal adhesion potential.

        Parameters:
            d: distance
            dhat_p: distance of largest adhesion force (:math:`\hat{d}_p`) where :math:`0 < \hat{d}_p < \hat{d}_a`
            dhat_a: adhesion activation distance (:math:`\hat{d}_a`)
            a2: adjustable parameter relating to the maximum derivative of a (:math:`a_2`)

        Returns:
            The normal adhesion potential.
        )ipc_Qu8mg5v7",
        "d"_a, "dhat_p"_a, "dhat_a"_a, "a2"_a);

    m.def(
        "normal_adhesion_potential_first_derivative",
        &normal_adhesion_potential_first_derivative,
        R"ipc_Qu8mg5v7(
        The first derivative of the normal adhesion potential wrt d.

        Parameters:
            d: distance
            dhat_p: distance of largest adhesion force (:math:`\hat{d}_p`) where :math:`0 < \hat{d}_p < \hat{d}_a`
            dhat_a: adhesion activation distance (:math:`\hat{d}_a`)
            a2: adjustable parameter relating to the maximum derivative of a (:math:`a_2`)

        Returns:
            The first derivative of the normal adhesion potential wrt d.
        )ipc_Qu8mg5v7",
        "d"_a, "dhat_p"_a, "dhat_a"_a, "a2"_a);

    m.def(
        "normal_adhesion_potential_second_derivative",
        &normal_adhesion_potential_second_derivative,
        R"ipc_Qu8mg5v7(
        The second derivative of the normal adhesion potential wrt d.

        Parameters:
            d: distance
            dhat_p: distance of largest adhesion force (:math:`\hat{d}_p`) where :math:`0 < \hat{d}_p < \hat{d}_a`
            dhat_a: adhesion activation distance (:math:`\hat{d}_a`)
            a2: adjustable parameter relating to the maximum derivative of a (:math:`a_2`)

        Returns:
            The second derivative of the normal adhesion potential wrt d.
        )ipc_Qu8mg5v7",
        "d"_a, "dhat_p"_a, "dhat_a"_a, "a2"_a);

    m.def(
        "max_normal_adhesion_force_magnitude",
        &max_normal_adhesion_force_magnitude,
        R"ipc_Qu8mg5v7(
        The maximum normal adhesion force magnitude.

        Parameters:
            dhat_p: distance of largest adhesion force (:math:`\hat{d}_p`) where :math:`0 < \hat{d}_p < \hat{d}_a`
            dhat_a: adhesion activation distance (:math:`\hat{d}_a`)
            a2: adjustable parameter relating to the maximum derivative of a (:math:`a_2`)

        Returns:
            The maximum normal adhesion force magnitude.
        )ipc_Qu8mg5v7",
        "dhat_p"_a, "dhat_a"_a, "a2"_a);

    m.def(
        "tangential_adhesion_f0", &tangential_adhesion_f0,
        R"ipc_Qu8mg5v7(
        The tangential adhesion mollifier function.

        Parameters:
            y: The tangential relative speed.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The tangential adhesion mollifier function at y.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_a"_a);

    m.def(
        "tangential_adhesion_f1", &tangential_adhesion_f1,
        R"ipc_Qu8mg5v7(
        The first derivative of the tangential adhesion mollifier function.

        Parameters:
            y: The tangential relative speed.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The first derivative of the tangential adhesion mollifier function at y.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_a"_a);

    m.def(
        "tangential_adhesion_f2", &tangential_adhesion_f2,
        R"ipc_Qu8mg5v7(
        The second derivative of the tangential adhesion mollifier function.

        Parameters:
            y: The tangential relative speed.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The second derivative of the tangential adhesion mollifier function at y.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_a"_a);

    m.def(
        "tangential_adhesion_f1_over_x", &tangential_adhesion_f1_over_x,
        R"ipc_Qu8mg5v7(
        The first derivative of the tangential adhesion mollifier function divided by y.

        Parameters:
            y: The tangential relative speed.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The first derivative of the tangential adhesion mollifier function divided by y.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_a"_a);

    m.def(
        "tangential_adhesion_f2_x_minus_f1_over_x3",
        &tangential_adhesion_f2_x_minus_f1_over_x3,
        R"ipc_Qu8mg5v7(
        The second derivative of the tangential adhesion mollifier function times y minus the first derivative all divided by y cubed.

        Parameters:
            y: The tangential relative speed.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The second derivative of the tangential adhesion mollifier function times y minus the first derivative all divided by y cubed.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_a"_a);

    m.def(
        "smooth_mu_a0", &smooth_mu_a0,
        R"ipc_Qu8mg5v7(
        Compute the value of the ∫ μ(y) a₁(y) dy, where a₁ is the first derivative of the smooth tangential adhesion mollifier.

        Note:
            The `a0`/`a1` are unrelated to the `a0`/`a1` in the normal adhesion.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static adhesion.
            mu_k: Coefficient of kinetic adhesion.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The value of the integral at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_a"_a);

    m.def(
        "smooth_mu_a1", &smooth_mu_a1,
        R"ipc_Qu8mg5v7(
        Compute the value of the μ(y) a₁(y), where a₁ is the first derivative of the smooth tangential adhesion mollifier.

        Note:
            The `a1` is unrelated to the `a1` in the normal adhesion.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static adhesion.
            mu_k: Coefficient of kinetic adhesion.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The value of the product at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_a"_a);

    m.def(
        "smooth_mu_a2", &smooth_mu_a2,
        R"ipc_Qu8mg5v7(
        Compute the value of d/dy (μ(y) a₁(y)), where a₁ is the first derivative of the smooth tangential adhesion mollifier.

        Note:
            The `a1`/`a2` are unrelated to the `a1`/`a2` in the normal adhesion.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static adhesion.
            mu_k: Coefficient of kinetic adhesion.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The value of the derivative at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_a"_a);

    m.def(
        "smooth_mu_a1_over_x", &smooth_mu_a1_over_x,
        R"ipc_Qu8mg5v7(
        Compute the value of the μ(y) a₁(y) / y, where a₁ is the first derivative of the smooth tangential adhesion mollifier.

        Notes:
            The `x` in the function name refers to the parameter `y`.
            The `a1` is unrelated to the `a1` in the normal adhesion.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static adhesion.
            mu_k: Coefficient of kinetic adhesion.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The value of the product at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_a"_a);

    m.def(
        "smooth_mu_a2_x_minus_mu_a1_over_x3",
        &smooth_mu_a2_x_minus_mu_a1_over_x3,
        R"ipc_Qu8mg5v7(
        Compute the value of the [(d/dy μ(y) a₁(y)) ⋅ y - μ(y) a₁(y)] / y³, where a₁ and a₂ are the first and second derivatives of the smooth tangential adhesion mollifier.

        Notes:
            The `x` in the function name refers to the parameter `y`.
            The `a1`/`a2` are unrelated to the `a1`/`a2` in the normal adhesion.

        Parameters:
            y: The tangential relative speed.
            mu_s: Coefficient of static adhesion.
            mu_k: Coefficient of kinetic adhesion.
            eps_a: Velocity threshold below which static adhesion force is applied.

        Returns:
            The value of the expression at y.
        )ipc_Qu8mg5v7",
        "y"_a, "mu_s"_a, "mu_k"_a, "eps_a"_a);
}
