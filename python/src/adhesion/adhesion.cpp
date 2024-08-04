#include <common.hpp>

#include <ipc/adhesion/adhesion.hpp>

namespace py = pybind11;
using namespace ipc;

void define_adhesion(py::module_& m)
{
    m.def(
        "normal_adhesion_potential", &normal_adhesion_potential,
        R"ipc_Qu8mg5v7(
        The normal adhesion potential.

        Parameters:
            d: distance
            dhat_p
            dhat_a
            a2

        Returns:

        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat_p"), py::arg("dhat_a"), py::arg("a2"));

    m.def(
        "normal_adhesion_potential_first_derivative",
        &normal_adhesion_potential_first_derivative,
        R"ipc_Qu8mg5v7(
        The first derivative of the normal adhesion potential wrt d.

        Parameters:
            d: distance
            dhat_p
            dhat_a
            a2

        Returns:

        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat_p"), py::arg("dhat_a"), py::arg("a2"));

    m.def(
        "normal_adhesion_potential_second_derivative",
        &normal_adhesion_potential_second_derivative,
        R"ipc_Qu8mg5v7(
        The second derivative of the normal adhesion potential wrt d.

        Parameters:
            d: distance
            dhat_p
            dhat_a
            a2

        Returns:

        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat_p"), py::arg("dhat_a"), py::arg("a2"));

    m.def(
        "f0_t", &f0_t,
        R"ipc_Qu8mg5v7(
        The tangent adhesion potential.

        Parameters:
            y
            eps_a

        Returns:

        )ipc_Qu8mg5v7",
        py::arg("y"), py::arg("eps_a"));

    m.def(
        "f1_t", &f1_t, "The tangent adhesion potential gradient.", py::arg("y"),
        py::arg("eps_a"));

    m.def("df1_t", &df1_t, py::arg("y"), py::arg("eps_a"));

    m.def("f1_t_over_x", &f1_t_over_x, py::arg("y"), py::arg("eps_a"));

    m.def(
        "df1_t_x_minus_f1_t_over_x3", &df1_t_x_minus_f1_t_over_x3, py::arg("y"),
        py::arg("eps_a"));
}
