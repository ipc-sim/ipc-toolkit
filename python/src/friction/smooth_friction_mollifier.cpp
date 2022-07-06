#include "../common.hpp"

#include <ipc/friction/smooth_friction_mollifier.hpp>

namespace py = pybind11;
using namespace ipc;

void define_smooth_friction_mollifier(py::module_& m)
{
    m.def("f0_SF", &f0_SF<double>, "", py::arg("x"), py::arg("epsv_times_h"));

    m.def(
        "f1_SF_over_x", &f1_SF_over_x<double>,
        "Derivative of f0_SF divided by x", py::arg("x"),
        py::arg("epsv_times_h"));

    m.def(
        "df1_x_minus_f1_over_x3", &df1_x_minus_f1_over_x3<double>,
        "$\frac{f_1'(x)x + f_1(x)}{x^3}$", py::arg("x"),
        py::arg("epsv_times_h"));
}
