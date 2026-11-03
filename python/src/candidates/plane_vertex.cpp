#include <ipc/candidates/plane_vertex.hpp>

#include <common.hpp>

#include "stencil_methods.hpp"

using namespace ipc;

void define_plane_vertex_candidate(py::module_& m)
{
    auto cls = py::class_<PlaneVertexCandidate>(m, "PlaneVertexCandidate");
    cls.def(
           py::init<Eigen::Hyperplane<double, 3>, index_t>(), "plane"_a,
           "vertex_id"_a)
        .def_readwrite(
            "plane", &PlaneVertexCandidate::plane, "Plane of the candidate")
        .def_readwrite(
            "vertex_id", &PlaneVertexCandidate::vertex_id, "ID of the vertex");

    define_stencil_methods<PlaneVertexCandidate>(cls);

    // The adapter the collision types derive from; it carries the
    // polymorphic stencil interface the POD candidate deliberately
    // lacks.
    py::class_<PlaneVertexStencil, PlaneVertexCandidate, CollisionStencil>(
        m, "PlaneVertexStencil");
}