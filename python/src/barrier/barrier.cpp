#include <common.hpp>

#include <ipc/barrier/barrier.hpp>

namespace py = pybind11;
using namespace ipc;

class PyBarrier : public Barrier {
public:
    // clang-format off
    using Barrier::Barrier; // Inherit constructors
    double operator()(const double d, const double dhat) const override { PYBIND11_OVERRIDE_PURE_NAME(double, Barrier, "__call__", operator(), d, dhat); }
    double first_derivative(const double d, const double dhat) const override { PYBIND11_OVERRIDE_PURE(double, Barrier, first_derivative, d, dhat); }
    double second_derivative(const double d, const double dhat) const override { PYBIND11_OVERRIDE_PURE(double, Barrier, second_derivative, d, dhat); }
    // clang-format on
};

void define_barrier(py::module_& m)
{
    py::class_<Barrier, PyBarrier, std::shared_ptr<Barrier>>(m, "Barrier")
        .def(py::init<>())
        .def("__call__", &Barrier::operator(), py::arg("d"), py::arg("dhat"))
        .def(
            "first_derivative", &Barrier::first_derivative, py::arg("d"),
            py::arg("dhat"))
        .def(
            "second_derivative", &Barrier::second_derivative, py::arg("d"),
            py::arg("dhat"));

    py::class_<ClampedLogBarrier, Barrier, std::shared_ptr<ClampedLogBarrier>>(
        m, "ClampedLogBarrier")
        .def(py::init());

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
        "barrier_first_derivative", &barrier_first_derivative,
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
        "barrier_second_derivative", &barrier_second_derivative,
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
