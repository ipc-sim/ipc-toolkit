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

    // New pairwise functions
    m.def(
        "f0_SF_pairwise", &f0_SF_pairwise,
        R"ipc_Qu8mg5v7(
        Smooth friction mollifier function for pairwise materials.

        .. math::

            f_0(s)= \mu \begin{cases}
            -\frac{s^3}{3\epsilon_v^2} + \frac{s^2}{\epsilon_v}
            + \frac{\epsilon_v}{3}, & |s| < \epsilon_v
            \newline
            s, & |s| \geq \epsilon_v
            \end{cases}

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            mu1: Friction coefficient for material 1.
            mu2: Friction coefficient for material 2.

        Returns:
            The value of the mollifier function at s, using the blended mu.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("mu1"), py::arg("mu2"));

    m.def(
        "f1_SF_over_x_pairwise", &f1_SF_over_x_pairwise,
        R"ipc_Qu8mg5v7(
        Compute the derivative of f0_SF_pairwise divided by s (:math:`\frac{f_0'(s)}{s}`).

        .. math::

            f_1(s) = \mu f_0'(s) = \begin{cases}
            -\frac{s^2}{\epsilon_v^2}+\frac{2 s}{\epsilon_v}, & |s| < \epsilon_v
            \newline 1, & |s| \geq \epsilon_v
            \end{cases}

        .. math::

            \frac{f_1(s)}{s} = \mu \begin{cases}
            -\frac{s}{\epsilon_v^2}+\frac{2}{\epsilon_v}, & |s| < \epsilon_v
            \newline \frac{1}{s}, & |s| \geq \epsilon_v
            \end{cases}

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            mu1: Friction coefficient for material 1.
            mu2: Friction coefficient for material 2.

        Returns:
            The value of the derivative of f0_SF divided by s, using blended mu.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("mu1"), py::arg("mu2"));

    m.def(
        "df1_x_minus_f1_over_x3_pairwise", &df1_x_minus_f1_over_x3_pairwise,
        R"ipc_Qu8mg5v7(
        The derivative of f1 times s minus f1 all divided by s cubed, using blended mu.

        .. math::

            \frac{f_1'(s) s - f_1(s)}{s^3} = \mu \begin{cases}
            -\frac{1}{s \epsilon_v^2}, & |s| < \epsilon_v \newline
            -\frac{1}{s^3}, & |s| \geq \epsilon_v
            \end{cases}

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            mu1: Friction coefficient for material 1.
            mu2: Friction coefficient for material 2.

        Returns:
            The derivative of f1 times s minus f1 all divided by s cubed, using blended mu.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("mu1"), py::arg("mu2"));

    // New pairwise transition functions
    m.def(
        "f0_SF_pairwise_transition", &f0_SF_pairwise_transition,
        R"ipc_Qu8mg5v7(
        Smooth friction mollifier function for pairwise materials with transition between static and kinetic friction.

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            static_mu: Static friction coefficient.
            kinetic_mu: Kinetic friction coefficient.

        Returns:
            The value of the mollifier function at s, using the appropriate mu (static or kinetic).
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("static_mu"), py::arg("kinetic_mu"));

    m.def(
        "f1_SF_over_x_pairwise_transition", &f1_SF_over_x_pairwise_transition,
        R"ipc_Qu8mg5v7(
        Compute the derivative of f0_SF_pairwise_transition divided by s (:math:`\frac{f_0'(s)}{s}`).

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            static_mu: Static friction coefficient.
            kinetic_mu: Kinetic friction coefficient.

        Returns:
            The value of the derivative of f0_SF divided by s, using the appropriate mu (static or kinetic).
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("static_mu"), py::arg("kinetic_mu"));

    m.def(
        "df1_x_minus_f1_over_x3_pairwise_transition", &df1_x_minus_f1_over_x3_pairwise_transition,
        R"ipc_Qu8mg5v7(
        The derivative of f1 times s minus f1 all divided by s cubed, using the appropriate mu (static or kinetic).

        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            static_mu: Static friction coefficient.
            kinetic_mu: Kinetic friction coefficient.

        Returns:
            The derivative of f1 times s minus f1 all divided by s cubed, using the appropriate mu (static or kinetic).
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("static_mu"), py::arg("kinetic_mu"));
}
