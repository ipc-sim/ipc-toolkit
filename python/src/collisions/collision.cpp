#include <common.hpp>

#include <ipc/collisions/collision.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collision(py::module_& m)
{
    py::class_<Collision<4>, CollisionStencil<4>>(m, "Collision")
        .def(
            "is_mollified", &Collision<4>::is_mollified,
            "Does the distance potentially have to be mollified?")
        .def(
            "mollifier_threshold", &Collision<4>::mollifier_threshold,
            R"ipc_Qu8mg5v7(
            Compute the mollifier threshold for the distance.

            Parameters:
                rest_positions: The stencil's rest vertex positions.

            Returns:
                The mollifier threshold.
            )ipc_Qu8mg5v7",
            py::arg("rest_positions"))
        .def(
            "mollifier",
            py::overload_cast<const VectorMax12d&>(
                &Collision<4>::mollifier, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the mollifier for the distance.

            Parameters:
                positions: The stencil's vertex positions.

            Returns:
                The mollifier value.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "mollifier",
            py::overload_cast<const VectorMax12d&, double>(
                &Collision<4>::mollifier, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the mollifier for the distance.

            Parameters:
                positions: The stencil's vertex positions.
                eps_x: The mollifier's threshold.

            Returns:
                The mollifier value.
            )ipc_Qu8mg5v7",
            py::arg("positions"), py::arg("eps_x"))
        .def(
            "mollifier_gradient",
            py::overload_cast<const VectorMax12d&>(
                &Collision<4>::mollifier_gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the mollifier for the distance wrt the positions.

            Parameters:
                positions: The stencil's vertex positions.

            Returns:
                The mollifier gradient.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "mollifier_gradient",
            py::overload_cast<const VectorMax12d&, double>(
                &Collision<4>::mollifier_gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the mollifier for the distance wrt the positions.

            Parameters:
                positions: The stencil's vertex positions.
                eps_x: The mollifier's threshold.

            Returns:
                The mollifier gradient.
            )ipc_Qu8mg5v7",
            py::arg("positions"), py::arg("eps_x"))
        .def(
            "mollifier_hessian",
            py::overload_cast<const VectorMax12d&>(
                &Collision<4>::mollifier_hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the Hessian of the mollifier for the distance wrt the positions.

            Parameters:
                positions: The stencil's vertex positions.

            Returns:
                The mollifier Hessian.
            )ipc_Qu8mg5v7",
            py::arg("positions"))
        .def(
            "mollifier_hessian",
            py::overload_cast<const VectorMax12d&, double>(
                &Collision<4>::mollifier_hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the Hessian of the mollifier for the distance wrt the positions.

            Parameters:
                positions: The stencil's vertex positions.
                eps_x: The mollifier's threshold.

            Returns:
                The mollifier Hessian.
            )ipc_Qu8mg5v7",
            py::arg("positions"), py::arg("eps_x"))
        .def(
            "mollifier_gradient_wrt_x", &Collision<4>::mollifier_gradient_wrt_x,
            R"ipc_Qu8mg5v7(
            Compute the gradient of the mollifier for the distance w.r.t. rest positions.

            Parameters:
                rest_positions: The stencil's rest vertex positions.
                positions: The stencil's vertex positions.

            Returns:
                The mollifier gradient w.r.t. rest positions.
            )ipc_Qu8mg5v7",
            py::arg("rest_positions"), py::arg("positions"))
        .def(
            "mollifier_gradient_jacobian_wrt_x",
            &Collision<4>::mollifier_gradient_jacobian_wrt_x,
            R"ipc_Qu8mg5v7(
            Compute the jacobian of the distance mollifier's gradient w.r.t. rest positions.

            Parameters:
                rest_positions: The stencil's rest vertex positions.
                positions: The stencil's vertex positions.

            Returns:
                The jacobian of the mollifier's gradient w.r.t. rest positions.
            )ipc_Qu8mg5v7",
            py::arg("rest_positions"), py::arg("positions"))
        .def_readwrite(
            "dmin", &Collision<4>::dmin, "The minimum separation distance.")
        .def_readwrite(
            "weight", &Collision<4>::weight,
            "The term's weight (e.g., collision area)")
        .def_property(
            "weight_gradient",
            [](const Collision<4>& self) -> Eigen::SparseMatrix<double> {
                return self.weight_gradient;
            },
            [](Collision<4>& self,
               const Eigen::SparseMatrix<double>& weight_gradient) {
                assert_is_sparse_vector(weight_gradient, "weight_gradient");
                self.weight_gradient = weight_gradient;
            },
            "The gradient of the term's weight wrt the rest positions.");
}
