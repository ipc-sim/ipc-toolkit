#include <common.hpp>

#include <ipc/collisions/smooth_point_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_smooth_point_edge(py::module_& m)
{
    m.def(
        "smooth_point_edge_potential_pointwise", &smooth_point_edge_potential_pointwise,
        R"ipc_Qu8mg5v7(
        Compute the new potential between a point and a point on an edge.

        Parameters:
            p, e0, e1, uv, dhat, alpha

        Returns:
            potential

        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"), py::arg("uv"),
        py::arg("dhat"), py::arg("alpha"));

    m.def(
        "smooth_point_edge_potential", &smooth_point_edge_potential,
        R"ipc_Qu8mg5v7(
        Compute the new potential between a point and an edge.

        Parameters:
            p, e0, e1, dhat, alpha

        Returns:
            potential

        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dhat"), py::arg("alpha"));
}
