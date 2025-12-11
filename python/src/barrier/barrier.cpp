#include <common.hpp>

#include <ipc/barrier/barrier.hpp>

using namespace ipc;

class PyBarrier : public Barrier {
public:
    // clang-format off
    using Barrier::Barrier; // Inherit constructors
    double operator()(const double d, const double dhat) const override { PYBIND11_OVERRIDE_PURE_NAME(double, Barrier, "__call__", operator(), d, dhat); }
    double first_derivative(const double d, const double dhat) const override { PYBIND11_OVERRIDE_PURE(double, Barrier, first_derivative, d, dhat); }
    double second_derivative(const double d, const double dhat) const override { PYBIND11_OVERRIDE_PURE(double, Barrier, second_derivative, d, dhat); }
    double units(const double dhat) const override { PYBIND11_OVERRIDE_PURE(double, Barrier, units, dhat); }
    // clang-format on
};

void define_barrier(py::module_& m)
{
    py::class_<Barrier, PyBarrier, std::shared_ptr<Barrier>>(m, "Barrier")
        .def(py::init<>())
        .def(
            "__call__", &Barrier::operator(), "d"_a, "dhat"_a,
            R"ipc_Qu8mg5v7(
            Evaluate the barrier function.

            Parameters:
                d: The distance.
                dhat: Activation distance of the barrier.

            Returns:
                The value of the barrier function at d.
            )ipc_Qu8mg5v7")
        .def(
            "first_derivative", &Barrier::first_derivative, "d"_a, "dhat"_a,
            R"ipc_Qu8mg5v7(
            Evaluate the first derivative of the barrier function wrt d.

            Parameters:
                d: The distance.
                dhat: Activation distance of the barrier.

            Returns:
                The value of the first derivative of the barrier function at d.
            )ipc_Qu8mg5v7")
        .def(
            "second_derivative", &Barrier::second_derivative, "d"_a, "dhat"_a,
            R"ipc_Qu8mg5v7(
            Evaluate the second derivative of the barrier function wrt d.

            Parameters:
                d: The distance.
                dhat: Activation distance of the barrier.

            Returns:
                The value of the second derivative of the barrier function at d.
            )ipc_Qu8mg5v7")
        .def(
            "units", &Barrier::units, "dhat"_a,
            R"ipc_Qu8mg5v7(
            Get the units of the barrier function.

            Parameters:
                dhat: Activation distance of the barrier.

            Returns:
                The units of the barrier function.
            )ipc_Qu8mg5v7");

    py::class_<ClampedLogBarrier, Barrier, std::shared_ptr<ClampedLogBarrier>>(
        m, "ClampedLogBarrier",
        R"ipc_Qu8mg5v7(
        Smoothly clamped log barrier functions from [Li et al. 2020].

        .. math::

            b(d) = -(d-\hat{d})^2\ln\left(\frac{d}{\hat{d}}\right)

        )ipc_Qu8mg5v7")
        .def(py::init());

    py::class_<
        NormalizedClampedLogBarrier, ClampedLogBarrier,
        std::shared_ptr<NormalizedClampedLogBarrier>>(
        m, "NormalizedClampedLogBarrier",
        R"ipc_Qu8mg5v7(
        Normalized barrier function from [Li et al. 2023].

        .. math::

            b(d) = -\left(\frac{d}{\hat{d}}-1\right)^2\ln\left(\frac{d}{\hat{d}}\right)

        )ipc_Qu8mg5v7")
        .def(py::init());

    py::class_<
        ClampedLogSqBarrier, Barrier, std::shared_ptr<ClampedLogSqBarrier>>(
        m, "ClampedLogSqBarrier",
        R"ipc_Qu8mg5v7(
        Clamped log barrier with a quadratic log term from [Huang et al. 2024].

        .. math::

            b(d) = (d-\hat{d})^2\ln^2\left(\frac{d}{\hat{d}}\right)

        )ipc_Qu8mg5v7")
        .def(py::init());

    py::class_<CubicBarrier, Barrier, std::shared_ptr<CubicBarrier>>(
        m, "CubicBarrier", R"ipc_Qu8mg5v7(
        Cubic barrier function from [Ando 2024].

        .. math::

            b(d) = -\frac{2}{3\hat{d}} (d - \hat{d})^3

        )ipc_Qu8mg5v7")
        .def(py::init());

    py::class_<TwoStageBarrier, Barrier, std::shared_ptr<TwoStageBarrier>>(
        m, "TwoStageBarrier",
        R"ipc_Qu8mg5v7(
        Two-stage activation barrier from [Chen et al. 2025].

        .. math::

            b(d) = \begin{cases}
                -\frac{\hat{d}^2}{4} \left(\ln\left(\frac{2d}{\hat{d}}\right) -
                \tfrac{1}{2}\right) & d < \frac{\hat{d}}{2}\\
                \tfrac{1}{2} (\hat{d} - d)^2 & d < \hat{d}\\
                0 & d \ge \hat{d}
            \end{cases}

        )ipc_Qu8mg5v7")
        .def(py::init());

    m.def(
        "barrier", &barrier, "d"_a, "dhat"_a,
        R"ipc_Qu8mg5v7(
        Function that grows to infinity as d approaches 0 from the right.

        .. math::

            b(d) = -(d-\hat{d})^2\ln\left(\frac{d}{\hat{d}}\right)

        Parameters:
            d: The distance.
            dhat: Activation distance of the barrier.

        Returns:
            The value of the barrier function at d.
        )ipc_Qu8mg5v7");

    m.def(
        "barrier_first_derivative", &barrier_first_derivative, "d"_a, "dhat"_a,
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
        )ipc_Qu8mg5v7");

    m.def(
        "barrier_second_derivative", &barrier_second_derivative, "d"_a,
        "dhat"_a,
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
        )ipc_Qu8mg5v7");
}
