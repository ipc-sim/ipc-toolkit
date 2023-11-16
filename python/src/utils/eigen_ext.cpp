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
        R"ipc_Qu8mg5v7(
        Matrix projection onto positive definite cone

        Parameters:
            A: Symmetric matrix to project

        Returns:
            Projected matrix
        )ipc_Qu8mg5v7",
        py::arg("A"), py::arg("eps") = 1e-8);

    m.def(
        "project_to_psd",
        &project_to_psd<
            double, Eigen::Dynamic, Eigen::Dynamic,
            Eigen::ColMajor | Eigen::AutoAlign, Eigen::Dynamic, Eigen::Dynamic>,
        R"ipc_Qu8mg5v7(
        Matrix projection onto positive semi-definite cone

        Parameters:
            A: Symmetric matrix to project

        Returns:
            Projected matrix
        )ipc_Qu8mg5v7",
        py::arg("A"));
}
