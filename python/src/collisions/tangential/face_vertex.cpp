#include <common.hpp>

#include <ipc/collisions/tangential/face_vertex.hpp>

using namespace ipc;

void define_face_vertex_tangential_collision(py::module_& m)
{
    py::class_<
        FaceVertexTangentialCollision, FaceVertexCandidate,
        TangentialCollision>(m, "FaceVertexTangentialCollision")
        .def(py::init<const FaceVertexNormalCollision&>(), "collision"_a)
        .def(
            py::init<
                const FaceVertexNormalCollision&, Eigen::ConstRef<VectorMax12d>,
                const NormalPotential&, const double>(),
            "collision"_a, "positions"_a, "normal_potential"_a,
            "normal_stiffness"_a);
}
