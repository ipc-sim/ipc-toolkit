#include "../common.hpp"

#include <ipc/broad_phase/voxel_size_heuristic.hpp>

namespace py = pybind11;
using namespace ipc;

void define_voxel_size_heuristic(py::module_& m)
{
    m.def(
        "suggest_good_voxel_size",
        py::overload_cast<
            const Eigen::MatrixXd&, const Eigen::MatrixXi&, double>(
            &suggest_good_voxel_size),
        "", py::arg("V"), py::arg("E"), py::arg("inflation_radius") = 0);

    m.def(
        "suggest_good_voxel_size",
        py::overload_cast<
            const Eigen::MatrixXd&, const Eigen::MatrixXd&,
            const Eigen::MatrixXi&, double>(&suggest_good_voxel_size),
        "", py::arg("V0"), py::arg("V1"), py::arg("E"),
        py::arg("inflation_radius") = 0);

    m.def(
        "mean_edge_length",
        [](const Eigen::MatrixXd& V0, const Eigen::MatrixXd& V1,
           const Eigen::MatrixXi& E) {
            double std_deviation;
            double r = mean_edge_length(V0, V1, E, std_deviation);
            return std::make_tuple(r, std_deviation);
        },
        "Compute the average edge length of a mesh.", py::arg("V0"),
        py::arg("V1"), py::arg("E"));

    m.def(
        "mean_displacement_length",
        [](const Eigen::MatrixXd& U) {
            double std_deviation;
            double r = mean_displacement_length(U, std_deviation);
            return std::make_tuple(r, std_deviation);
        },
        "Compute the average displacement length.", py::arg("U"));

    m.def(
        "median_edge_length", &median_edge_length,
        "Compute the median edge length of a mesh.", py::arg("V0"),
        py::arg("V1"), py::arg("E"));

    m.def(
        "median_displacement_length", &median_displacement_length,
        "Compute the median displacement length.", py::arg("U"));
}
