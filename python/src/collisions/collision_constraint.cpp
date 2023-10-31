#include <common.hpp>

#include <ipc/collisions/collision_constraint.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision_constraint(py::module_& m)
{
    py::class_<CollisionConstraint, CollisionStencil>(m, "CollisionConstraint")
        // .def(
        //     py::init<const double, const Eigen::SparseVector<double>&>(),
        //     py::arg("weight"), py::arg("weight_gradient"))
        .def(
            "compute_potential", &CollisionConstraint::compute_potential,
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &CollisionConstraint::compute_potential_gradient,
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &CollisionConstraint::compute_potential_hessian,
            py::arg("vertices"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"), py::arg("project_hessian_to_psd"))
        .def(
            "compute_shape_derivative",
            [](const CollisionConstraint& self,
               const Eigen::MatrixXd& rest_positions,
               const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges,
               const Eigen::MatrixXi& faces, const double dhat) {
                std::vector<Eigen::Triplet<double>> triplets;
                self.compute_shape_derivative(
                    rest_positions, vertices, edges, faces, dhat, triplets);
                return triplets;
            },
            "Compute the derivative of the potential gradient wrt the shape.",
            py::arg("rest_positions"), py::arg("vertices"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"))
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
