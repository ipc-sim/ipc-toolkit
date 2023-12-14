#include <common.hpp>

#include <ipc/collisions/collision_constraint.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_constraint(py::module_& m)
{
    py::class_<CollisionConstraint, CollisionStencil>(m, "CollisionConstraint")
        .def_readwrite("dmin", &CollisionConstraint::dmin)
        .def_readwrite("weight", &CollisionConstraint::weight)
        .def_property(
            "weight_gradient",
            [](const CollisionConstraint& self) -> Eigen::SparseMatrix<double> {
                return self.weight_gradient;
            },
            [](CollisionConstraint& self,
               const Eigen::SparseMatrix<double>& weight_gradient) {
                assert_is_sparse_vector(weight_gradient, "weight_gradient");
                self.weight_gradient = weight_gradient;
            });
}
