#include <common.hpp>

#include <ipc/friction/collisions/friction_collision.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_collision(py::module_& m)
{
    py::class_<FrictionCollision, CollisionStencil>(m, "FrictionCollision")
        .def_readwrite(
            "normal_force_magnitude",
            &FrictionCollision::normal_force_magnitude,
            "Collision force magnitude")
        .def_readwrite("mu", &FrictionCollision::mu, "Coefficient of friction")
        .def_readwrite("weight", &FrictionCollision::weight, "Weight")
        .def_property(
            "weight_gradient",
            [](const FrictionCollision& self) -> Eigen::SparseMatrix<double> {
                return self.weight_gradient;
            },
            [](FrictionCollision& self,
               const Eigen::SparseMatrix<double>& weight_gradient) {
                assert_is_sparse_vector(weight_gradient, "weight_gradient");
                self.weight_gradient = weight_gradient;
            },
            "Gradient of weight with respect to all DOF")
        .def_readwrite(
            "closest_point", &FrictionCollision::closest_point,
            "Barycentric coordinates of the closest point(s)")
        .def_readwrite(
            "tangent_basis", &FrictionCollision::tangent_basis,
            "Tangent basis of the collision (max size 3×2)");
}
