#include <common.hpp>

#include <ipc/collisions/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_constraint(py::module_& m)
{
    py::class_<FaceVertexConstraint, FaceVertexCandidate, CollisionConstraint>(
        m, "FaceVertexConstraint")
        .def(
            py::init<long, long>(), "", py::arg("face_id"),
            py::arg("vertex_id"))
        .def(py::init<FaceVertexCandidate>(), "", py::arg("fv_candidate"));
}
