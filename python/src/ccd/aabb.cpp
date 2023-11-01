#include <common.hpp>

#include <ipc/ccd/aabb.hpp>

namespace py = pybind11;
using namespace ipc;

void define_ccd_aabb(py::module_& m)
{
    m.def(
        "point_edge_aabb_cd", &point_edge_aabb_cd, py::arg("p"), py::arg("e0"),
        py::arg("e1"), py::arg("dist"));

    m.def(
        "edge_edge_aabb_cd", &edge_edge_aabb_cd, py::arg("ea0"), py::arg("ea1"),
        py::arg("eb0"), py::arg("eb1"), py::arg("dist"));

    m.def(
        "point_triangle_aabb_cd", &point_triangle_aabb_cd, py::arg("p"),
        py::arg("t0"), py::arg("t1"), py::arg("t2"), py::arg("dist"));

    m.def(
        "edge_triangle_aabb_cd", &edge_triangle_aabb_cd, py::arg("e0"),
        py::arg("e1"), py::arg("t0"), py::arg("t1"), py::arg("t2"),
        py::arg("dist"));

    m.def(
        "point_edge_aabb_ccd", &point_edge_aabb_ccd, py::arg("p_t0"),
        py::arg("e0_t0"), py::arg("e1_t0"), py::arg("p_t1"), py::arg("e0_t1"),
        py::arg("e1_t1"), py::arg("dist"));

    m.def(
        "edge_edge_aabb_ccd", &edge_edge_aabb_ccd, py::arg("ea0_t0"),
        py::arg("ea1_t0"), py::arg("eb0_t0"), py::arg("eb1_t0"),
        py::arg("ea0_t1"), py::arg("ea1_t1"), py::arg("eb0_t1"),
        py::arg("eb1_t1"), py::arg("dist"));

    m.def(
        "point_triangle_aabb_ccd", &point_triangle_aabb_ccd, py::arg("p_t0"),
        py::arg("t0_t0"), py::arg("t1_t0"), py::arg("t2_t0"), py::arg("p_t1"),
        py::arg("t0_t1"), py::arg("t1_t1"), py::arg("t2_t1"), py::arg("dist"));
}
