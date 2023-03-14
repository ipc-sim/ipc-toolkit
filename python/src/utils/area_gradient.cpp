#include <common.hpp>

#include <ipc/utils/area_gradient.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace py = pybind11;
using namespace ipc;

void define_area_gradient(py::module_& m)
{
    m.def(
        "edge_length_gradient",
        [](const VectorMax3d& e0, const VectorMax3d& e1) {
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            VectorMax6d grad;
            edge_length_gradient(e0, e1, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of an edge's length.

        Parameters:
            e0: The first vertex of the edge.
            e1: The second vertex of the edge.

        Returns:
            The gradient of the edge's length wrt e0, and e1.
        )ipc_Qu8mg5v7",
        py::arg("e0"), py::arg("e1"));

    m.def(
        "triangle_area_gradient",
        [](const Eigen::Vector3d& t0, const Eigen::Vector3d& t1,
           const Eigen::Vector3d& t2) {
            Vector<double, 9> grad;
            triangle_area_gradient(t0, t1, t2, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the area of a triangle.

        Parameters:
            t0: The first vertex of the triangle.
            t1: The second vertex of the triangle.
            t2: The third vertex of the triangle.

        Returns:
            The gradient of the triangle's area t0, t1, and t2.
        )ipc_Qu8mg5v7",
        py::arg("t0"), py::arg("t1"), py::arg("t2"));
}
