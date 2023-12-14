#include <common.hpp>

#include <ipc/friction/constraints/friction_constraint.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_constraint(py::module_& m)
{
    py::class_<FrictionConstraint, CollisionStencil>(m, "FrictionConstraint")
        .def_readwrite(
            "normal_force_magnitude",
            &FrictionConstraint::normal_force_magnitude,
            "Contact force magnitude")
        .def_readwrite("mu", &FrictionConstraint::mu, "Coefficient of friction")
        .def_readwrite("weight", &FrictionConstraint::weight, "Weight")
        .def_property(
            "weight_gradient",
            [](const FrictionConstraint& self) -> Eigen::SparseMatrix<double> {
                return self.weight_gradient;
            },
            [](FrictionConstraint& self,
               const Eigen::SparseMatrix<double>& weight_gradient) {
                assert_is_sparse_vector(weight_gradient, "weight_gradient");
                self.weight_gradient = weight_gradient;
            },
            "Gradient of weight with respect to all DOF")
        .def_readwrite(
            "closest_point", &FrictionConstraint::closest_point,
            "Barycentric coordinates of the closest point(s)")
        .def_readwrite(
            "tangent_basis", &FrictionConstraint::tangent_basis,
            "Tangent basis of the contact (max size 3×2)");
}
