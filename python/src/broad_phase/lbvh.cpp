#include <common.hpp>

#include <ipc/broad_phase/lbvh.hpp>

using namespace ipc;

void define_lbvh(py::module_& m)
{
    py::class_<LBVH::Node>(m, "LBVH_Node")
        .def_readonly("aabb_min", &LBVH::Node::aabb_min)
        .def_readonly("aabb_max", &LBVH::Node::aabb_max)
        .def_readonly("left", &LBVH::Node::left)
        .def_readonly("right", &LBVH::Node::right)
        .def_property_readonly("is_leaf", &LBVH::Node::is_leaf)
        .def_property_readonly("is_inner", &LBVH::Node::is_inner)
        .def_property_readonly("is_valid", &LBVH::Node::is_valid);

    py::class_<LBVH, BroadPhase, std::shared_ptr<LBVH>>(m, "LBVH")
        .def(py::init())
        .def_property_readonly("vertex_nodes", &LBVH::vertex_nodes)
        .def_property_readonly("edge_nodes", &LBVH::edge_nodes)
        .def_property_readonly("face_nodes", &LBVH::face_nodes);
}
