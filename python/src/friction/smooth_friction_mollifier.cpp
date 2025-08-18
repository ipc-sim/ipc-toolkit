#include <common.hpp>

#include <ipc/friction/smooth_friction_mollifier.hpp>

using namespace ipc;

void define_smooth_friction_mollifier(py::module_& m)
{
    m.def(
        "smooth_friction_f0", &smooth_friction_f0,
        R"ipc_Qu8mg5v7(
        Smooth friction mollifier function.

        .. math::

            f_0(y)= \begin{cases}
            -\frac{y^3}{3\epsilon_v^2} + \frac{y^2}{\epsilon_v}
            + \frac{\epsilon_v}{3}, & |y| < \epsilon_v
            \newline
            y, & |y| \geq \epsilon_v
            \end{cases}

        Parameters:
            y: The tangential relative speed.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the mollifier function at y.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_v"_a);

    m.def(
        "smooth_friction_f1", &smooth_friction_f1,
        R"ipc_Qu8mg5v7(
        The first derivative of the smooth friction mollifier.

        .. math::

            f_1(y) = f_0'(y) = \begin{cases}
            -\frac{y^2}{\epsilon_v^2}+\frac{2 y}{\epsilon_v}, & |y| < \epsilon_v
            \newline 1, & |y| \geq \epsilon_v
            \end{cases}

        Parameters:
            y: The tangential relative speed.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the derivative of the smooth friction mollifier at y.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_v"_a);

    m.def(
        "smooth_friction_f2", &smooth_friction_f2,
        R"ipc_Qu8mg5v7(
        The second derivative of the smooth friction mollifier.

        .. math::

            f_2(y) = f_0''(y) = \begin{cases}
            -\frac{2 y}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |y| < \epsilon_v
            \newline 0, & |y| \geq \epsilon_v
            \end{cases}

        Parameters:
            y: The tangential relative speed.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the second derivative of the smooth friction mollifier at y.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_v"_a);

    m.def(
        "smooth_friction_f1_over_x", &smooth_friction_f1_over_x,
        R"ipc_Qu8mg5v7(
        Compute the derivative of the smooth friction mollifier divided by y (:math:`\frac{f_0'(y)}{y}`).

        .. math::

            \frac{f_1(y)}{y} = \begin{cases}
            -\frac{y}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |y| < \epsilon_v
            \newline \frac{1}{y}, & |y| \geq \epsilon_v
            \end{cases}

        Parameters:
            y: The tangential relative speed.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The value of the derivative of smooth_friction_f0 divided by y.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_v"_a);

    m.def(
        "smooth_friction_f2_x_minus_f1_over_x3",
        &smooth_friction_f2_x_minus_f1_over_x3,
        R"ipc_Qu8mg5v7(
        The derivative of f1 times y minus f1 all divided by y cubed.

        .. math::

            \frac{f_1'(y) y - f_1(y)}{y^3} = \begin{cases}
            -\frac{1}{y \epsilon_v^2}, & |y| < \epsilon_v \newline
            -\frac{1}{y^3}, & |y| \geq \epsilon_v
            \end{cases}

        Parameters:
            y: The tangential relative speed.
            eps_v: Velocity threshold below which static friction force is applied.

        Returns:
            The derivative of f1 times y minus f1 all divided by y cubed.
        )ipc_Qu8mg5v7",
        "y"_a, "eps_v"_a);
}
