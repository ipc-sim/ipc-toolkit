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
            dhat_p: distance of largest adhesion force (\f(\hat{d}_p\f)) (\f(0 < \hat{d}_p < \hat{d}_a\f))
            dhat_a: adhesion activation distance (\f(\hat{d}_a\f))
            a2: adjustable parameter relating to the maximum derivative of a (\f(a_2\f))

        Returns:
            The normal adhesion potential.
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat_p"), py::arg("dhat_a"), py::arg("a2"));

    m.def(
        "normal_adhesion_potential_first_derivative",
        &normal_adhesion_potential_first_derivative,
        R"ipc_Qu8mg5v7(
        The first derivative of the normal adhesion potential wrt d.

        Parameters:
            d: distance
            dhat_p: distance of largest adhesion force (\f(\hat{d}_p\f)) (\f(0 < \hat{d}_p < \hat{d}_a\f))
            dhat_a: adhesion activation distance (\f(\hat{d}_a\f))
            a2: adjustable parameter relating to the maximum derivative of a (\f(a_2\f))

        Returns:
            The first derivative of the normal adhesion potential wrt d.
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat_p"), py::arg("dhat_a"), py::arg("a2"));

    m.def(
        "normal_adhesion_potential_second_derivative",
        &normal_adhesion_potential_second_derivative,
        R"ipc_Qu8mg5v7(
        The second derivative of the normal adhesion potential wrt d.

        Parameters:
            d: distance
            dhat_p: distance of largest adhesion force (\f(\hat{d}_p\f)) (\f(0 < \hat{d}_p < \hat{d}_a\f))
            dhat_a: adhesion activation distance (\f(\hat{d}_a\f))
            a2: adjustable parameter relating to the maximum derivative of a (\f(a_2\f))

        Returns:
            The second derivative of the normal adhesion potential wrt d.
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat_p"), py::arg("dhat_a"), py::arg("a2"));

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
        py::arg("y"), py::arg("eps_a"));

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
        py::arg("y"), py::arg("eps_a"));

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
        py::arg("y"), py::arg("eps_a"));

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
        py::arg("y"), py::arg("eps_a"));

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
        py::arg("y"), py::arg("eps_a"));
}
