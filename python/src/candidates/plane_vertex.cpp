#include <common.hpp>

#include <ipc/candidates/plane_vertex.hpp>

using namespace ipc;

void define_plane_vertex_candidate(py::module_& m)
{
    py::class_<PlaneVertexCandidate, CollisionStencil>(
        m, "PlaneVertexCandidate")
        .def(
            py::init<
                Eigen::ConstRef<VectorMax3d>, Eigen::ConstRef<VectorMax3d>,
                index_t>(),
            "plane_origin"_a, "plane_normal"_a, "vertex_id"_a)
        .def_readwrite(
            "plane_origin", &PlaneVertexCandidate::plane_origin,
            "Origin of the plane")
        .def_readwrite(
            "plane_normal", &PlaneVertexCandidate::plane_normal,
            "Normal of the plane")
        .def_readwrite(
            "vertex_id", &PlaneVertexCandidate::vertex_id, "ID of the vertex");
}
