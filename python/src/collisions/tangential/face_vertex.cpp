#include <common.hpp>

#include <ipc/collisions/tangential/face_vertex.hpp>

namespace py = pybind11;
using namespace ipc;

void define_face_vertex_tangential_collision(py::module_& m)
{
    py::class_<
        FaceVertexTangentialCollision, FaceVertexCandidate,
        TangentialCollision>(m, "FaceVertexTangentialCollision")
        .def(py::init<const FaceVertexNormalCollision&>(), py::arg("collision"))
        .def(
            py::init<
                const FaceVertexNormalCollision&, Eigen::ConstRef<VectorMax12d>,
                const NormalPotential&, const double>(),
            py::arg("collision"), py::arg("positions"),
            py::arg("normal_potential"), py::arg("normal_stiffness"));
}
