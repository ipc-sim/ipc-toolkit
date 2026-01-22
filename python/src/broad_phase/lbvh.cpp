#include <common.hpp>

#include <ipc/broad_phase/lbvh.hpp>

using namespace ipc;

void define_lbvh(py::module_& m)
{
    py::class_<LBVH::Node>(m, "LBVH_Node")
        .def_readonly("aabb_min", &LBVH::Node::aabb_min)
        .def_readonly("aabb_max", &LBVH::Node::aabb_max)
        .def_property_readonly(
            "left",
            [](const LBVH::Node& node) {
                if (node.is_inner()) {
                    return node.left;
                } else {
                    throw py::attribute_error(
                        "Leaf nodes do not have 'left' child.");
                }
            })
        .def_property_readonly(
            "right",
            [](const LBVH::Node& node) {
                if (node.is_inner()) {
                    return node.right;
                } else {
                    throw py::attribute_error(
                        "Leaf nodes do not have 'right' child.");
                }
            })
        .def_property_readonly(
            "primitive_id",
            [](const LBVH::Node& node) {
                if (node.is_leaf()) {
                    return node.primitive_id;
                } else {
                    throw py::attribute_error(
                        "Inner nodes do not have a 'primitive_id'");
                }
            })
        .def_property_readonly("is_leaf", &LBVH::Node::is_leaf)
        .def_property_readonly("is_inner", &LBVH::Node::is_inner)
        .def_property_readonly("is_valid", &LBVH::Node::is_valid)
        .def("intersects", &LBVH::Node::intersects)
        .def("__str__", [](const LBVH::Node& node) {
            return fmt::format(
                "LBVH::Node(aabb_min={}, aabb_max={}, {})",
                node.aabb_min.transpose(), node.aabb_max.transpose(),
                node.is_inner()
                    ? fmt::format("left={}, right={}", node.left, node.right)
                    : fmt::format("primitive_id={}", node.primitive_id));
        });

    py::class_<LBVH, BroadPhase, std::shared_ptr<LBVH>>(m, "LBVH")
        .def(py::init())
        .def_property_readonly("vertex_nodes", &LBVH::vertex_nodes)
        .def_property_readonly("edge_nodes", &LBVH::edge_nodes)
        .def_property_readonly("face_nodes", &LBVH::face_nodes);
}
