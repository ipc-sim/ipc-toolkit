#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>

#include <ipc/utils/eigen_ext.hpp>

namespace py = pybind11;
using namespace ipc;

void define_eigen_ext_members(py::module_& m)
{
    m.def(
        "project_to_pd",
        &project_to_pd<
            double, Eigen::Dynamic, Eigen::Dynamic,
            Eigen::ColMajor | Eigen::AutoAlign, Eigen::Dynamic, Eigen::Dynamic>,
        "Matrix Projection onto Positive Definite Cone", py::arg("A"),
        py::arg("eps") = 1e-8);

    m.def(
        "project_to_psd",
        &project_to_psd<
            double, Eigen::Dynamic, Eigen::Dynamic,
            Eigen::ColMajor | Eigen::AutoAlign, Eigen::Dynamic, Eigen::Dynamic>,
        "Matrix Projection onto Positive Semi-Definite Cone", py::arg("A"));
}
