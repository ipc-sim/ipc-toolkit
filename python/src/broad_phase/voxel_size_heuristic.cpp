#include <common.hpp>

#include <ipc/broad_phase/voxel_size_heuristic.hpp>

using namespace ipc;

void define_voxel_size_heuristic(py::module_& m)
{
    m.def(
        "suggest_good_voxel_size",
        py::overload_cast<
            Eigen::ConstRef<Eigen::MatrixXd>, Eigen::ConstRef<Eigen::MatrixXi>,
            const double>(&suggest_good_voxel_size),
        "vertices"_a, "edges"_a, "inflation_radius"_a = 0);

    m.def(
        "suggest_good_voxel_size",
        py::overload_cast<
            Eigen::ConstRef<Eigen::MatrixXd>, Eigen::ConstRef<Eigen::MatrixXd>,
            Eigen::ConstRef<Eigen::MatrixXi>, const double>(
            &suggest_good_voxel_size),
        "vertices_t0"_a, "vertices_t1"_a, "edges"_a, "inflation_radius"_a = 0);

    m.def(
        "mean_edge_length",
        [](Eigen::ConstRef<Eigen::MatrixXd> vertices_t0,
           Eigen::ConstRef<Eigen::MatrixXd> vertices_t1,
           Eigen::ConstRef<Eigen::MatrixXi> edges) {
            double std_deviation;
            double r = mean_edge_length(
                vertices_t0, vertices_t1, edges, std_deviation);
            return std::make_tuple(r, std_deviation);
        },
        "Compute the average edge length of a mesh.", "vertices_t0"_a,
        "vertices_t1"_a, "edges"_a);

    m.def(
        "mean_displacement_length",
        [](Eigen::ConstRef<Eigen::MatrixXd> displacements) {
            double std_deviation;
            double r = mean_displacement_length(displacements, std_deviation);
            return std::make_tuple(r, std_deviation);
        },
        "Compute the average displacement length.", "displacements"_a);

    m.def(
        "median_edge_length", &median_edge_length,
        "Compute the median edge length of a mesh.", "vertices_t0"_a,
        "vertices_t1"_a, "edges"_a);

    m.def(
        "median_displacement_length", &median_displacement_length,
        "Compute the median displacement length.", "displacements"_a);

    m.def(
        "max_edge_length", &max_edge_length,
        "Compute the maximum edge length of a mesh.", "vertices_t0"_a,
        "vertices_t1"_a, "edges"_a);

    m.def(
        "max_displacement_length", &max_displacement_length,
        "Compute the maximum displacement length.", "displacements"_a);
}
