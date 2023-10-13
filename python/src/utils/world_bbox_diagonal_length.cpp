#include <common.hpp>

#include <ipc/utils/world_bbox_diagonal_length.hpp>

namespace py = pybind11;
using namespace ipc;

void define_world_bbox_diagonal_length(py::module_& m)
{
    m.def(
        "world_bbox_diagonal_length", &world_bbox_diagonal_length,
        R"ipc_Qu8mg5v7(
        Compute the diagonal length of the world bounding box.

        Parameters:
            vertices: Vertex positions

        Returns:
            The diagonal length of the world bounding box.
        )ipc_Qu8mg5v7",
        py::arg("vertices"));
}
