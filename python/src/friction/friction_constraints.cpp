#include "../common.hpp"

#include <ipc/friction/friction_constraints.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_constraints(py::module_& m)
{
    py::class_<FrictionConstraints>(m, "FrictionConstraints")
        .def(py::init(), "")
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const CollisionConstraints&, double, double, double>(
                &FrictionConstraints::build),
            "", py::arg("mesh"), py::arg("vertices"),
            py::arg("contact_constraints"), py::arg("dhat"),
            py::arg("barrier_stiffness"), py::arg("mu"))
        .def(
            "build",
            [](FrictionConstraints& self, const CollisionMesh& mesh,
               const Eigen::MatrixXd& vertices,
               const CollisionConstraints& contact_constraints, double dhat,
               double barrier_stiffness, const Eigen::VectorXd& mus) {
                self.build(
                    mesh, vertices, contact_constraints, dhat,
                    barrier_stiffness, mus);
            },
            "", py::arg("mesh"), py::arg("vertices"),
            py::arg("contact_constraints"), py::arg("dhat"),
            py::arg("barrier_stiffness"), py::arg("mus"))
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const CollisionConstraints&, double, double,
                const Eigen::VectorXd&,
                const std::function<double(double, double)>&>(
                &FrictionConstraints::build),
            "", py::arg("mesh"), py::arg("vertices"),
            py::arg("contact_constraints"), py::arg("dhat"),
            py::arg("barrier_stiffness"), py::arg("mus"), py::arg("blend_mu"))
        .def(
            "compute_potential",
            &FrictionConstraints::compute_potential<double>,
            R"ipc_Qu8mg5v7(
            Compute the friction dissapative potential from the given velocity.

            Parameters:
                mesh: The collision mesh.
                velocity: Current vertex velocity (rowwise).

            Returns:
                The friction dissapative potential.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("velocity"), py::arg("epsv_times_h"))
        .def(
            "compute_potential_gradient",
            &FrictionConstraints::compute_potential_gradient,
            R"ipc_Qu8mg5v7(
            Compute the gradient of the friction dissapative potential wrt the velocity.

            Parameters:
                mesh: The collision mesh.
                velocity: Current vertex velocity (rowwise).

            Returns:
                The gradient of the friction dissapative potential wrt the velocity.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("velocity"), py::arg("epsv_times_h"))
        .def(
            "compute_potential_hessian",
            &FrictionConstraints::compute_potential_hessian,
            R"ipc_Qu8mg5v7(
            Compute the Hessian of the friction dissapative potential wrt the velocity.

            Parameters:
                mesh: The collision mesh.
                velocity: Current vertex velocity (rowwise).
                project_hessian_to_psd: If true, project the Hessian to be positive semi-definite.

            Returns:
                The Hessian of the friction dissapative potential wrt the velocity.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("velocity"), py::arg("epsv_times_h"),
            py::arg("project_hessian_to_psd") = true)
        .def(
            "compute_force",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const Eigen::MatrixXd&, const double,
                const double, const double, const double, const bool>(
                &FrictionConstraints::compute_force, py::const_),
            "", py::arg("mesh"), py::arg("X"), py::arg("Ut"), py::arg("U"),
            py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("dmin") = 0,
            py::arg("no_mu") = false)
        .def(
            "compute_force",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const double, const double,
                const double, const double, const bool>(
                &FrictionConstraints::compute_force, py::const_),
            "", py::arg("mesh"), py::arg("X"), py::arg("U"), py::arg("dhat"),
            py::arg("barrier_stiffness"), py::arg("epsv_times_h"),
            py::arg("dmin") = 0, py::arg("no_mu") = false)
        .def(
            "compute_force_jacobian",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const Eigen::MatrixXd&, const double,
                const double, const double, const FrictionConstraint::DiffWRT,
                const double>(
                &FrictionConstraints::compute_force_jacobian, py::const_),
            "", py::arg("mesh"), py::arg("X"), py::arg("Ut"), py::arg("U"),
            py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("wrt"), py::arg("dmin") = 0)
        .def(
            "compute_force_jacobian",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const double, const double,
                const double, const FrictionConstraint::DiffWRT, const double>(
                &FrictionConstraints::compute_force_jacobian, py::const_),
            "", py::arg("mesh"), py::arg("X"), py::arg("U"), py::arg("dhat"),
            py::arg("barrier_stiffness"), py::arg("epsv_times_h"),
            py::arg("wrt"), py::arg("dmin") = 0)
        .def(
            "__len__", &FrictionConstraints::size,
            "Get the number of friction constraints.")
        .def(
            "empty", &FrictionConstraints::empty,
            "Get if the friction constraints are empty.")
        .def(
            "clear", &FrictionConstraints::clear,
            "Clear the friction constraints.")
        .def(
            "__getitem__",
            [](FrictionConstraints& self, size_t idx) -> FrictionConstraint& {
                return self[idx];
            },
            py::return_value_policy::reference,
            R"ipc_Qu8mg5v7(
            Get a reference to constriant idx.

            Parameters:
                idx: The index of the constraint.

            Returns:
                A reference to the constraint.
            )ipc_Qu8mg5v7",
            py::arg("idx"))
        .def_readwrite(
            "vv_constraints", &FrictionConstraints::vv_constraints, "")
        .def_readwrite(
            "ev_constraints", &FrictionConstraints::ev_constraints, "")
        .def_readwrite(
            "ee_constraints", &FrictionConstraints::ee_constraints, "")
        .def_readwrite(
            "fv_constraints", &FrictionConstraints::fv_constraints, "");
}
