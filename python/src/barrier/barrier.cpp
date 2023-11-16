#include <common.hpp>

#include <ipc/barrier/barrier.hpp>

namespace py = pybind11;
using namespace ipc;

void define_barrier(py::module_& m)
{
    m.def(
        "barrier", &barrier,
        R"ipc_Qu8mg5v7(
        Function that grows to infinity as d approaches 0 from the right.

        .. math::

            b(d) = -(d-\hat{d})^2\ln\left(\frac{d}{\hat{d}}\right)

        Parameters:
            d: The distance.
            dhat: Activation distance of the barrier.

        Returns:
            The value of the barrier function at d.
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));

    m.def(
        "barrier_gradient", &barrier_gradient,
        R"ipc_Qu8mg5v7(
        Derivative of the barrier function.

        .. math::

            b'(d) = (\hat{d}-d) \left( 2\ln\left( \frac{d}{\hat{d}} \right) -
            \frac{\hat{d}}{d} + 1\right)

        Parameters:
            d: The distance.
            dhat: Activation distance of the barrier.

        Returns:
            The derivative of the barrier wrt d.
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));

    m.def(
        "barrier_hessian", &barrier_hessian,
        R"ipc_Qu8mg5v7(
        Second derivative of the barrier function.

        .. math::

            b''(d) = \left( \frac{\hat{d}}{d} + 2 \right) \frac{\hat{d}}{d} -
            2\ln\left( \frac{d}{\hat{d}} \right) - 3

        Parameters:
            d: The distance.
            dhat: Activation distance of the barrier.

        Returns:
            The second derivative of the barrier wrt d.
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));
}
