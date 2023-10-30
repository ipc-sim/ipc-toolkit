#include <common.hpp>

#include <ipc/broad_phase/voxel_size_heuristic.hpp>

namespace py = pybind11;
using namespace ipc;

void define_voxel_size_heuristic(py::module_& m)
{
    m.def(
        "suggest_good_voxel_size",
        py::overload_cast<
            const Eigen::MatrixXd&, const Eigen::MatrixXi&, const double>(
            &suggest_good_voxel_size),
        py::arg("vertices"), py::arg("edges"), py::arg("inflation_radius") = 0);

    m.def(
        "suggest_good_voxel_size",
        py::overload_cast<
            const Eigen::MatrixXd&, const Eigen::MatrixXd&,
            const Eigen::MatrixXi&, const double>(&suggest_good_voxel_size),
        py::arg("vertices_t0"), py::arg("vertices_t1"), py::arg("edges"),
        py::arg("inflation_radius") = 0);

    m.def(
        "mean_edge_length",
        [](const Eigen::MatrixXd& vertices_t0,
           const Eigen::MatrixXd& vertices_t1, const Eigen::MatrixXi& edges) {
            double std_deviation;
            double r = mean_edge_length(
                vertices_t0, vertices_t1, edges, std_deviation);
            return std::make_tuple(r, std_deviation);
        },
        "Compute the average edge length of a mesh.", py::arg("vertices_t0"),
        py::arg("vertices_t1"), py::arg("edges"));

    m.def(
        "mean_displacement_length",
        [](const Eigen::MatrixXd& displacements) {
            double std_deviation;
            double r = mean_displacement_length(displacements, std_deviation);
            return std::make_tuple(r, std_deviation);
        },
        "Compute the average displacement length.", py::arg("displacements"));

    m.def(
        "median_edge_length", &median_edge_length,
        "Compute the median edge length of a mesh.", py::arg("vertices_t0"),
        py::arg("vertices_t1"), py::arg("edges"));

    m.def(
        "median_displacement_length", &median_displacement_length,
        "Compute the median displacement length.", py::arg("displacements"));

    m.def(
        "max_edge_length", &max_edge_length,
        "Compute the maximum edge length of a mesh.", py::arg("vertices_t0"),
        py::arg("vertices_t1"), py::arg("edges"));

    m.def(
        "max_displacement_length", &max_displacement_length,
        "Compute the maximum displacement length.", py::arg("displacements"));
}
