#include "../common.hpp"

#include <ipc/friction/friction.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction(py::module_& m)
{
    m.def(
        "construct_friction_constraint_set",
        [](const CollisionMesh& mesh, const Eigen::MatrixXd& V,
           const Constraints& contact_constraint_set, double dhat,
           double barrier_stiffness, double mu) {
            FrictionConstraints friction_constraint_set;
            construct_friction_constraint_set(
                mesh, V, contact_constraint_set, dhat, barrier_stiffness, mu,
                friction_constraint_set);
            return friction_constraint_set;
        },
        "", py::arg("mesh"), py::arg("V"), py::arg("contact_constraint_set"),
        py::arg("dhat"), py::arg("barrier_stiffness"), py::arg("mu"));

    m.def(
        "construct_friction_constraint_set",
        [](const CollisionMesh& mesh, const Eigen::MatrixXd& V,
           const Constraints& contact_constraint_set, double dhat,
           double barrier_stiffness, const Eigen::VectorXd& mus) {
            FrictionConstraints friction_constraint_set;
            construct_friction_constraint_set(
                mesh, V, contact_constraint_set, dhat, barrier_stiffness, mus,
                friction_constraint_set);
            return friction_constraint_set;
        },
        "", py::arg("mesh"), py::arg("V"), py::arg("contact_constraint_set"),
        py::arg("dhat"), py::arg("barrier_stiffness"), py::arg("mus"));

    m.def(
        "construct_friction_constraint_set",
        [](const CollisionMesh& mesh, const Eigen::MatrixXd& V,
           const Constraints& contact_constraint_set, double dhat,
           double barrier_stiffness, const Eigen::VectorXd& mus,
           const std::function<double(double, double)>& blend_mu) {
            FrictionConstraints friction_constraint_set;
            construct_friction_constraint_set(
                mesh, V, contact_constraint_set, dhat, barrier_stiffness, mus,
                blend_mu, friction_constraint_set);
            return friction_constraint_set;
        },
        "", py::arg("mesh"), py::arg("V"), py::arg("contact_constraint_set"),
        py::arg("dhat"), py::arg("barrier_stiffness"), py::arg("mus"),
        py::arg("blend_mu"));

    m.def(
        "compute_friction_potential", &compute_friction_potential<double>,
        R"ipc_Qu8mg5v7(
        Compute the friction potential between to positions.

        Parameters:
            V0: Vertex positions at start of time-step (rowwise)
            V1: Current vertex positions (rowwise)
            E:  Edge vertex indicies
            F:  Face vertex indicies (empty in 2D)
            friction_constraint_set
            epsv_times_h
        )ipc_Qu8mg5v7",
        py::arg("mesh"), py::arg("V0"), py::arg("V1"),
        py::arg("friction_constraint_set"), py::arg("epsv_times_h"));

    m.def(
        "compute_friction_potential_gradient",
        &compute_friction_potential_gradient, "", py::arg("mesh"),
        py::arg("V0"), py::arg("V1"), py::arg("friction_constraint_set"),
        py::arg("epsv_times_h"));

    m.def(
        "compute_friction_potential_hessian",
        &compute_friction_potential_hessian, "", py::arg("mesh"), py::arg("V0"),
        py::arg("V1"), py::arg("friction_constraint_set"),
        py::arg("epsv_times_h"), py::arg("project_hessian_to_psd") = true);

    m.def(
        "compute_friction_force",
        py::overload_cast<
            const CollisionMesh&, const Eigen::MatrixXd&,
            const Eigen::MatrixXd&, const Eigen::MatrixXd&,
            const FrictionConstraints&, const double, const double,
            const double, const double, const bool>(&compute_friction_force),
        "", py::arg("mesh"), py::arg("X"), py::arg("Ut"), py::arg("U"),
        py::arg("friction_constraint_set"), py::arg("dhat"),
        py::arg("barrier_stiffness"), py::arg("epsv_times_h"),
        py::arg("dmin") = 0, py::arg("no_mu") = false);

    m.def(
        "compute_friction_force",
        py::overload_cast<
            const CollisionMesh&, const Eigen::MatrixXd&,
            const Eigen::MatrixXd&, const FrictionConstraints&, const double,
            const double, const double, const double, const bool>(
            &compute_friction_force),
        "", py::arg("mesh"), py::arg("X"), py::arg("U"),
        py::arg("friction_constraint_set"), py::arg("dhat"),
        py::arg("barrier_stiffness"), py::arg("epsv_times_h"),
        py::arg("dmin") = 0, py::arg("no_mu") = false);

    m.def(
        "compute_friction_force_jacobian",
        py::overload_cast<
            const CollisionMesh&, const Eigen::MatrixXd&,
            const Eigen::MatrixXd&, const Eigen::MatrixXd&,
            const FrictionConstraints&, const double, const double,
            const double, const FrictionConstraint::DiffWRT, const double>(
            &compute_friction_force_jacobian),
        "", py::arg("mesh"), py::arg("X"), py::arg("Ut"), py::arg("U"),
        py::arg("friction_constraint_set"), py::arg("dhat"),
        py::arg("barrier_stiffness"), py::arg("epsv_times_h"), py::arg("wrt"),
        py::arg("dmin") = 0);

    m.def(
        "compute_friction_force_jacobian",
        py::overload_cast<
            const CollisionMesh&, const Eigen::MatrixXd&,
            const Eigen::MatrixXd&, const FrictionConstraints&, const double,
            const double, const double, const FrictionConstraint::DiffWRT,
            const double>(&compute_friction_force_jacobian),
        "", py::arg("mesh"), py::arg("X"), py::arg("U"),
        py::arg("friction_constraint_set"), py::arg("dhat"),
        py::arg("barrier_stiffness"), py::arg("epsv_times_h"), py::arg("wrt"),
        py::arg("dmin") = 0);
}
