#include <common.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace py = pybind11;
using namespace ipc;

void define_eigen_ext(py::module_& m)
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
