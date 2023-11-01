#include <common.hpp>

#include <ipc/friction/smooth_friction_mollifier.hpp>

namespace py = pybind11;
using namespace ipc;

void define_smooth_friction_mollifier(py::module_& m)
{
    m.def(
        "f0_SF", &f0_SF,
        R"ipc_Qu8mg5v7(
        Smooth friction mollifier function.

        .. math::

            f_0(s)= \begin{cases}
            -\frac{s^3}{3\epsilon_v^2} + \frac{s^2}{\epsilon_v}
            + \frac{\epsilon_v}{3}, & |s| < \epsilon_v
            \newline
            s, & |s| \geq \epsilon_v
            \end{cases}

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.

        Returns:
            The value of the mollifier function at s.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"));

    m.def(
        "f1_SF_over_x", &f1_SF_over_x,
        R"ipc_Qu8mg5v7(
        Compute the derivative of f0_SF divided by s (:math:`\frac{f_0'(s)}{s}`).

        .. math::

            f_1(s) = f_0'(s) = \begin{cases}
            -\frac{s^2}{\epsilon_v^2}+\frac{2 s}{\epsilon_v}, & |s| < \epsilon_v
            \newline 1, & |s| \geq \epsilon_v
            \end{cases}

        .. math::

            \frac{f_1(s)}{s} = \begin{cases}
            -\frac{s}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |s| < \epsilon_v
            \newline \frac{1}{s}, & |s| \geq \epsilon_v
            \end{cases}

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.

        Returns:
            The value of the derivative of f0_SF divided by s.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"));

    m.def(
        "df1_x_minus_f1_over_x3", &df1_x_minus_f1_over_x3,
        R"ipc_Qu8mg5v7(
        The derivative of f1 times s minus f1 all divided by s cubed.

        .. math::

            \frac{f_1'(s) s - f_1(s)}{s^3} = \begin{cases}
            -\frac{1}{s \epsilon_v^2}, & |s| < \epsilon_v \newline
            -\frac{1}{s^3}, & |s| \geq \epsilon_v
            \end{cases}

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.

        Returns:
            The derivative of f1 times s minus f1 all divided by s cubed.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"));
}
