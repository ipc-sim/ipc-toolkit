#include <common.hpp>

#include <ipc/utils/eigen_ext.hpp>

namespace py = pybind11;
using namespace ipc;

void define_eigen_ext(py::module_& m)
{
    m.def(
        "project_to_pd",
        // clang-format off
        &project_to_pd<
            double, Eigen::Dynamic, Eigen::Dynamic,
            Eigen::ColMajor | Eigen::AutoAlign, Eigen::Dynamic, Eigen::Dynamic>,
        // clang-format on
        R"ipc_Qu8mg5v7(
        Matrix projection onto positive definite cone

        Parameters:
            A: Symmetric matrix to project

        Returns:
            Projected matrix
        )ipc_Qu8mg5v7",
        py::arg("A"), py::arg("eps") = 1e-8);

    py::enum_<PSDProjectionMethod>(
        m, "PSDProjectionMethod",
        "Enumeration of implemented PSD projection methods.")
        .value("NONE", PSDProjectionMethod::NONE, "No PSD projection")
        .value(
            "CLAMP", PSDProjectionMethod::CLAMP,
            "Clamp negative eigenvalues to zero")
        .value("ABS", PSDProjectionMethod::ABS, "Flip negative eigenvalues")
        .export_values();

    m.def(
        "project_to_psd",
        // clang-format off
        &project_to_psd<
            double, Eigen::Dynamic, Eigen::Dynamic,
            Eigen::ColMajor | Eigen::AutoAlign, Eigen::Dynamic, Eigen::Dynamic>,
        // clang-format on
        R"ipc_Qu8mg5v7(
        Matrix projection onto positive semi-definite cone

        Parameters:
            A: Symmetric matrix to project
            method: PSD projection method

        Returns:
            Projected matrix
        )ipc_Qu8mg5v7",
        py::arg("A"), py::arg("method") = PSDProjectionMethod::CLAMP);
}
