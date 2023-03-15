#include <common.hpp>

#include <ipc/friction/smooth_friction_mollifier.hpp>

namespace py = pybind11;
using namespace ipc;

void define_smooth_friction_mollifier(py::module_& m)
{
    m.def(
        "f0_SF", &f0_SF<double>,
        R"ipc_Qu8mg5v7(
        Smooth friction mollifier function.

        .. math::

            f_0(y)= \begin{cases}
            -\frac{x^3}{3\epsilon_v^2 h^2} + \frac{x^2}{\epsilon_vh}
            + \frac{\epsilon_v h}{3}, & |x| \in\left[0, \epsilon_v h\right)
            \newline
            x, & |x| \geq \epsilon_v h
            \end{cases}

        Parameters:
            x: The tangential relative speed.
            epsv_times_h: Mollifier parameter :math:`\epsilon_v h\f`.

        Returns:
            The value of the mollifier function at x.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("epsv_times_h"));

    m.def(
        "f1_SF_over_x", &f1_SF_over_x<double>,
        R"ipc_Qu8mg5v7(
        Compute the derivative of f0_SF divided by x (:math:`\frac{f_0'(x)}{x}\f`).

        .. math::

            f_1(x) = f_0'(x) = \begin{cases}
            -\frac{x^2}{\epsilon_v^2 h^2}+\frac{2 x}{\epsilon_v h}, & |x|
            \in\left[0, \epsilon_v h \right) \newline
            1, & |x| \geq h \epsilon_v
            \end{cases}

        .. math::

            \frac{f_1(x)}{x} = \begin{cases}
            -\frac{x}{\epsilon_v^2 h^2}+\frac{2}{\epsilon_v h}, & |x|
            \in\left[0, \epsilon_v h \right) \newline
            \frac{1}{x}, & |x| \geq h \epsilon_v
            \end{cases}

        Parameters:
            x: The tangential relative speed.
            epsv_times_h: Mollifier parameter :math:`\epsilon_v h\f`.

        Returns:
            The value of the derivative of f0_SF divided by x.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("epsv_times_h"));

    m.def(
        "df1_x_minus_f1_over_x3", &df1_x_minus_f1_over_x3<double>,
        R"ipc_Qu8mg5v7(
        The derivative of f1 times x minus f1 all divided by x cubed.

        .. math::

            \frac{f_1'(x) x - f_1(x)}{x^3} = \begin{cases}
            -\frac{1}{x \epsilon_v^2 h^2}, & |x|
            \in\left[0, \epsilon_v h \right) \newline
            -\frac{1}{x^3}, & |x| \geq h \epsilon_v
            \end{cases}

        Parameters:
            x: The tangential relative speed.
            epsv_times_h: Mollifier parameter :math:`\epsilon_v h\f`.

        Returns:
            The derivative of f1 times x minus f1 all divided by x cubed.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("epsv_times_h"));
}
