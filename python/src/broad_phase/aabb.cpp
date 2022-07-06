#include <pybind11/pybind11.h>

#include <ipc/broad_phase/aabb.hpp>

namespace py = pybind11;
using namespace ipc;

void define_aabb_members(py::module_& m)
{
    py::class_<AABB>(m, "AABB")
        .def(py::init(), "")
        .def(
            py::init<const ArrayMax3d&, const ArrayMax3d&>(), "",
            py::arg("min"), py::arg("max"))
        .def(
            py::init<const AABB&, const AABB&>(), "", py::arg("aabb1"),
            py::arg("aabb2"))
        .def(
            py::init<const AABB&, const AABB&, const AABB&>(), "",
            py::arg("aabb1"), py::arg("aabb2"), py::arg("aabb3"))
        .def_static(
            "from_point",
            py::overload_cast<const VectorMax3d&, double>(&AABB::from_point),
            "Compute a AABB for a static point.", py::arg("p"),
            py::arg("inflation_radius") = 0)
        .def_static(
            "from_point",
            py::overload_cast<const VectorMax3d&, const VectorMax3d&, double>(
                &AABB::from_point),
            "Compute a AABB for a moving point (i.e. temporal edge).",
            py::arg("p_t0"), py::arg("p_t1"), py::arg("inflation_radius") = 0)
        .def("intersects", &AABB::intersects, "", py::arg("other"))
        .def_readwrite("min", &AABB::min, "")
        .def_readwrite("max", &AABB::max, "")
        .def_readwrite("vertex_ids", &AABB::vertex_ids, "");

    m.def(
        "build_vertex_boxes",
        [](const Eigen::MatrixXd& V, double inflation_radius = 0) {
            std::vector<AABB> vertex_boxes;
            build_vertex_boxes(V, vertex_boxes, inflation_radius);
            return vertex_boxes;
        },
        "", py::arg("V"), py::arg("inflation_radius") = 0);

    m.def(
        "build_vertex_boxes",
        [](const Eigen::MatrixXd& V0, const Eigen::MatrixXd& V1,
           double inflation_radius = 0) {
            std::vector<AABB> vertex_boxes;
            build_vertex_boxes(V0, V1, vertex_boxes, inflation_radius);
            return vertex_boxes;
        },
        "", py::arg("V0"), py::arg("V1"), py::arg("inflation_radius") = 0);

    m.def(
        "build_edge_boxes",
        [](const std::vector<AABB>& vertex_boxes, const Eigen::MatrixXi& E) {
            std::vector<AABB> edge_boxes;
            build_edge_boxes(vertex_boxes, E, edge_boxes);
            return edge_boxes;
        },
        "", py::arg("vertex_boxes"), py::arg("E"));

    m.def(
        "build_face_boxes",
        [](const std::vector<AABB>& vertex_boxes, const Eigen::MatrixXi& F) {
            std::vector<AABB> face_boxes;
            build_face_boxes(vertex_boxes, F, face_boxes);
            return face_boxes;
        },
        "", py::arg("vertex_boxes"), py::arg("F"));
}
