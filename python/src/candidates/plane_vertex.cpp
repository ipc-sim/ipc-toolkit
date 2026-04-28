#include <ipc/candidates/plane_vertex.hpp>

#include <common.hpp>

using namespace ipc;

void define_plane_vertex_candidate(py::module_& m)
{
    py::class_<PlaneVertexCandidate, CollisionStencil>(
        m, "PlaneVertexCandidate")
        .def(
            py::init<Eigen::Hyperplane<double, 3>, index_t>(), "plane"_a,
            "vertex_id"_a)
        .def_readwrite(
            "plane", &PlaneVertexCandidate::plane, "Plane of the candidate")
        .def_readwrite(
            "vertex_id", &PlaneVertexCandidate::vertex_id, "ID of the vertex");
}