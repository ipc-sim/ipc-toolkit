#include <common.hpp>

#include <ipc/collisions/collision.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision(py::module_& m)
{
    py::class_<Collision, CollisionStencil>(m, "Collision")
        .def_readwrite("dmin", &Collision::dmin)
        .def_readwrite("weight", &Collision::weight)
        .def_property(
            "weight_gradient",
            [](const Collision& self) -> Eigen::SparseMatrix<double> {
                return self.weight_gradient;
            },
            [](Collision& self,
               const Eigen::SparseMatrix<double>& weight_gradient) {
                assert_is_sparse_vector(weight_gradient, "weight_gradient");
                self.weight_gradient = weight_gradient;
            });
}
