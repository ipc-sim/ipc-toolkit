#include <common.hpp>

#include <ipc/geometry/angle.hpp>

using namespace ipc;

void define_angle(py::module_& m)
{
    m.def(
        "dihedral_angle", &dihedral_angle,
        R"ipc_Qu8mg5v7(
            Compute the bending angle between two triangles sharing an edge.
                x0---x2
                 | \ |
                x1---x3

            Parameters
            ----------
            x0 : Eigen::Vector3d
                The first vertex of the edge.
            x1 : Eigen::Vector3d
                The second vertex of the edge.
            x2 : Eigen::Vector3d
                The opposite vertex of the first triangle.
            x3 : Eigen::Vector3d
                The opposite vertex of the second triangle.

            Returns
            -------
            double
                The bending angle between the two triangles.
        )ipc_Qu8mg5v7",
        py::arg("x0"), py::arg("x1"), py::arg("x2"), py::arg("x3"));

    m.def(
        "dihedral_angle_gradient", &dihedral_angle_gradient,
        R"ipc_Qu8mg5v7(
            Compute the Jacobian of the bending angle between two triangles sharing an edge.
                x0---x2
                 | \ |
                x1---x3

            Parameters
            ----------
            x0 : Eigen::Vector3d
                The first vertex of the edge.
            x1 : Eigen::Vector3d
                The second vertex of the edge.
            x2 : Eigen::Vector3d
                The opposite vertex of the first triangle.
            x3 : Eigen::Vector3d
                The opposite vertex of the second triangle.

            Returns
            -------
            Eigen::Vector<double, 12>
                The Jacobian matrix of the bending angle with respect to the input vertices.
        )ipc_Qu8mg5v7",
        py::arg("x0"), py::arg("x1"), py::arg("x2"), py::arg("x3"));
}
