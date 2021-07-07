#include <pybind11/pybind11.h>

#include <ipc/barrier/barrier.hpp>

namespace py = pybind11;
using namespace ipc;

void define_barrier_functions(py::module_& m)
{
    m.def(
        "barrier", &barrier<double>,
        R"ipc_Qu8mg5v7(
        Function that grows to infinity as d approaches 0 from the right.

        Parameters
        ----------
        d : The distance
        dhat : Activation distance of the barrier

        Returns
        -------
        The value of the barrier function at d.

        Notes
        -----
        .. math:: b(d) = -(d-\hat{d})^2\ln(\frac{d}{\hat{d}})
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));

    m.def(
        "barrier_gradient", &barrier_gradient,
        R"ipc_Qu8mg5v7(
        Derivative of the barrier function.

        Parameters
        ----------
        d : The distance
        dhat : Activation distance of the barrier

        Returns
        -------
        The derivative of the barrier wrt d.

        Notes
        -----
        .. math:: b(d) = (\hat{d}-d)(2\ln(\frac{d}{\hat{d}}) - \frac{\hat{d}}{d}) + 1)
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));

    m.def(
        "barrier_hessian", &barrier_hessian,
        R"ipc_Qu8mg5v7(
        Second derivative of the barrier function.

        Parameters
        ----------
        d : The distance
        dhat : Activation distance of the barrier

        Returns
        -------
        The second derivative of the barrier wrt d.

        Notes
        -----
        .. math:: b(d) = (\frac{\hat{d}}{d} + 2)\frac{\hat{d}}{d} - 2\ln(\frac{d}{\hat{d}}) - 3
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));
}
