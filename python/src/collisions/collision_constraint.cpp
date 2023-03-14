#include <common.hpp>

#include <ipc/collisions/collision_constraint.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_constraint(py::module_& m)
{
    py::class_<CollisionConstraint, CollisionStencil>(m, "CollisionConstraint")
        .def(
            "compute_potential", &CollisionConstraint::compute_potential, "",
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &CollisionConstraint::compute_potential_gradient, "",
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &CollisionConstraint::compute_potential_hessian, "",
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"), py::arg("project_hessian_to_psd"))
        .def_readwrite(
            "minimum_distance", &CollisionConstraint::minimum_distance, "")
        .def_readwrite("weight", &CollisionConstraint::weight, "")
        .def_property(
            "weight_gradient",
            [](const CollisionConstraint& self) -> Eigen::SparseMatrix<double> {
                return self.weight_gradient;
            },
            [](CollisionConstraint& self,
               const Eigen::SparseMatrix<double>& weight_gradient) {
                assert_is_sparse_vector(weight_gradient, "weight_gradient");
                self.weight_gradient = weight_gradient;
            },
            "");
}
