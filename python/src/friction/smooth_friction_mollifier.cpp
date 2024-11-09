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
        ...
        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            static_mu: Static friction coefficient (optional).
            kinetic_mu: Kinetic friction coefficient (optional).
            blend_type: Type of blending to use (optional).

        Returns:
            The value of the mollifier function at s.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("static_mu") = std::nullopt,
        py::arg("kinetic_mu") = std::nullopt, py::arg("blend_type") = std::nullopt);

    m.def(
        "f1_SF_over_x", &f1_SF_over_x,
        R"ipc_Qu8mg5v7(
        Compute the derivative of f0_SF divided by s.
        ...
        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            static_mu: Static friction coefficient (optional).
            kinetic_mu: Kinetic friction coefficient (optional).
            blend_type: Type of blending to use (optional).

        Returns:
            The value of the derivative of f0_SF divided by s.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("static_mu") = std::nullopt,
        py::arg("kinetic_mu") = std::nullopt, py::arg("blend_type") = std::nullopt);

    m.def(
        "df1_x_minus_f1_over_x3", &df1_x_minus_f1_over_x3,
        R"ipc_Qu8mg5v7(
        The derivative of f1 times s minus f1 all divided by s cubed.
        ...
        Parameters:
            s: The tangential relative speed.
            epsv: Mollifier parameter :math:`\epsilon_v`.
            static_mu: Static friction coefficient (optional).
            kinetic_mu: Kinetic friction coefficient (optional).
            blend_type: Type of blending to use (optional).

        Returns:
            The derivative of f1 times s minus f1 all divided by s cubed.
        )ipc_Qu8mg5v7",
        py::arg("s"), py::arg("epsv"), py::arg("static_mu") = std::nullopt,
        py::arg("kinetic_mu") = std::nullopt, py::arg("blend_type") = std::nullopt);
}
