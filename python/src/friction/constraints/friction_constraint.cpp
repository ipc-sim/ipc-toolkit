#include <common.hpp>

#include <ipc/friction/constraints/friction_constraint.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_constraint(py::module_& m)
{
    py::class_<FrictionConstraint, CollisionStencil> friction_constraint(
        m, "FrictionConstraint");

    py::enum_<FrictionConstraint::DiffWRT>(friction_constraint, "DiffWRT")
        .value("X", FrictionConstraint::DiffWRT::X)
        .value("Ut", FrictionConstraint::DiffWRT::Ut)
        .value("U", FrictionConstraint::DiffWRT::U)
        .export_values();

    friction_constraint
        .def(
            "compute_potential_gradient",
            &FrictionConstraint::compute_potential_gradient,
            R"ipc_Qu8mg5v7(
            Compute the friction dissapative potential gradient wrt velocities.

            Parameters:
                velocities: Velocities of the vertices (rowwise)
                edges: Edges of the mesh
                faces: Faces of the mesh
                epsv_times_h: $\epsilon_vh$

            Returns:
                Gradient of the friction dissapative potential wrt velocities
            )ipc_Qu8mg5v7",
            py::arg("velocities"), py::arg("edges"), py::arg("faces"),
            py::arg("epsv_times_h"))
        .def(
            "compute_potential_hessian",
            &FrictionConstraint::compute_potential_hessian,
            R"ipc_Qu8mg5v7(
            Compute the friction dissapative potential hessian wrt velocities.

            Parameters:
                velocities: Velocities of the vertices (rowwise)
                edges: Edges of the mesh
                faces: Faces of the mesh
                epsv_times_h: $\epsilon_vh$
                project_hessian_to_psd: Project the hessian to PSD

            Returns:
                Hessian of the friction dissapative potential wrt velocities
            )ipc_Qu8mg5v7",
            py::arg("velocities"), py::arg("edges"), py::arg("faces"),
            py::arg("epsv_times_h"), py::arg("project_hessian_to_psd"))
        .def(
            "compute_force",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double, const double, const double, const bool>(
                &FrictionConstraint::compute_force, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force.

            Parameters:
                X: Rest positions of the vertices (rowwise)
                U: Current displacements of the vertices (rowwise)
                edges: Collision mesh edges
                faces: Collision mesh faces
                dhat: Barrier activation distance
                barrier_stiffness: Barrier stiffness
                epsv_times_h: $\epsilon_vh$
                dmin: Minimum distance
                no_mu: Whether to not multiply by mu

            Returns:
                Friction force
            )ipc_Qu8mg5v7",
            py::arg("X"), py::arg("U"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("dmin") = 0,
            py::arg("no_mu") = false)
        .def(
            "compute_force",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, const double, const double,
                const double, const double, const bool>(
                &FrictionConstraint::compute_force, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force.

            Parameters:
                X: Rest positions of the vertices (rowwise)
                Ut: Previous displacements of the vertices (rowwise)
                U: Current displacements of the vertices (rowwise)
                edges: Collision mesh edges
                faces: Collision mesh faces
                dhat: Barrier activation distance
                barrier_stiffness: Barrier stiffness
                epsv_times_h: $\epsilon_vh$
                dmin: Minimum distance
                no_mu: Whether to not multiply by mu

            Returns:
                Friction force
            )ipc_Qu8mg5v7",
            py::arg("X"), py::arg("Ut"), py::arg("U"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("dmin") = 0,
            py::arg("no_mu") = false)
        .def(
            "compute_force_jacobian",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXi&, const Eigen::MatrixXi&, const double,
                const double, const double, FrictionConstraint::DiffWRT,
                const double>(
                &FrictionConstraint::compute_force_jacobian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force Jacobian.

            Parameters:
                X: Rest positions of the vertices (rowwise)
                U: Current displacements of the vertices (rowwise)
                edges: Collision mesh edges
                faces: Collision mesh faces
                dhat: Barrier activation distance
                barrier_stiffness: Barrier stiffness
                epsv_times_h: $\epsilon_vh$
                wrt: Variable to differentiate the friction force with respect to.
                dmin: Minimum distance

            Returns:
                Friction force Jacobian
            )ipc_Qu8mg5v7",
            py::arg("X"), py::arg("U"), py::arg("edges"), py::arg("faces"),
            py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("wrt"), py::arg("dmin") = 0)
        .def(
            "compute_force_jacobian",
            py::overload_cast<
                const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                const Eigen::MatrixXd&, const Eigen::MatrixXi&,
                const Eigen::MatrixXi&, const double, const double,
                const double, FrictionConstraint::DiffWRT, const double>(
                &FrictionConstraint::compute_force_jacobian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force Jacobian.

            Parameters:
                X: Rest positions of the vertices (rowwise)
                Ut: Previous displacements of the vertices (rowwise)
                U: Current displacements of the vertices (rowwise)
                edges: Collision mesh edges
                faces: Collision mesh faces
                dhat: Barrier activation distance
                barrier_stiffness: Barrier stiffness
                epsv_times_h: $\epsilon_vh$
                wrt: Variable to differentiate the friction force with respect to.
                dmin: Minimum distance

            Returns:
                Friction force Jacobian
            )ipc_Qu8mg5v7",
            py::arg("X"), py::arg("Ut"), py::arg("U"), py::arg("edges"),
            py::arg("faces"), py::arg("dhat"), py::arg("barrier_stiffness"),
            py::arg("epsv_times_h"), py::arg("wrt"), py::arg("dmin") = 0)
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
