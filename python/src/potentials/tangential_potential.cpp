#include <common.hpp>
#include <potentials/potential.hpp>

#include <ipc/potentials/tangential_potential.hpp>

using namespace ipc;

void define_tangential_potential(py::module_& m)
{
    py::class_<TangentialPotential> tangential_potential(
        m, "TangentialPotential");

    define_potential_methods<TangentialCollisions>(tangential_potential);

    py::enum_<TangentialPotential::DiffWRT>(tangential_potential, "DiffWRT")
        .value(
            "REST_POSITIONS", TangentialPotential::DiffWRT::REST_POSITIONS,
            "Differentiate w.r.t. rest positions")
        .value(
            "LAGGED_DISPLACEMENTS",
            TangentialPotential::DiffWRT::LAGGED_DISPLACEMENTS,
            "Differentiate w.r.t. lagged displacements")
        .value(
            "VELOCITIES", TangentialPotential::DiffWRT::VELOCITIES,
            "Differentiate w.r.t. current velocities")
        .export_values();

    tangential_potential
        .def(
            "force",
            py::overload_cast<
                const TangentialCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXd>, const NormalPotential&,
                const double, const double, const bool>(
                &TangentialPotential::force, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force from the given velocities.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                rest_positions: Rest positions of the vertices (rowwise).
                lagged_displacements: Previous displacements of the vertices (rowwise).
                velocities: Current displacements of the vertices (rowwise).
                normal_potential: Normal potential (used for normal force magnitude).
                normal_stiffness: Normal stiffness (used for normal force magnitude).
                dmin: Minimum distance (used for normal force magnitude).
                no_mu: whether to not multiply by mu

            Returns:
                The friction force.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "rest_positions"_a,
            "lagged_displacements"_a, "velocities"_a, "normal_potential"_a,
            "normal_stiffness"_a, "dmin"_a = 0, "no_mu"_a = false)
        .def(
            "force_jacobian",
            py::overload_cast<
                const TangentialCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXd>,
                Eigen::ConstRef<Eigen::MatrixXd>, const NormalPotential&,
                const double, const TangentialPotential::DiffWRT, const double>(
                &TangentialPotential::force_jacobian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the Jacobian of the friction force wrt the velocities.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                rest_positions: Rest positions of the vertices (rowwise).
                lagged_displacements: Previous displacements of the vertices (rowwise).
                velocities: Current displacements of the vertices (rowwise).
                normal_potential: Normal potential (used for normal force magnitude).
                normal_stiffness: Normal stiffness (used for normal force magnitude).
                wrt: The variable to take the derivative with respect to.
                dmin: Minimum distance (used for normal force magnitude).

            Returns:
                The Jacobian of the friction force wrt the velocities.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "rest_positions"_a,
            "lagged_displacements"_a, "velocities"_a, "normal_potential"_a,
            "normal_stiffness"_a, "wrt"_a, "dmin"_a = 0)
        .def(
            "force",
            py::overload_cast<
                const TangentialCollision&, Eigen::ConstRef<VectorMax12d>,
                Eigen::ConstRef<VectorMax12d>, Eigen::ConstRef<VectorMax12d>,
                const NormalPotential&, const double, const double, const bool>(
                &TangentialPotential::force, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force.

            Parameters:
                collision: The collision
                rest_positions: Rest positions of the vertices (rowwise).
                lagged_displacements: Previous displacements of the vertices (rowwise).
                velocities: Current displacements of the vertices (rowwise).
                normal_potential: Normal potential (used for normal force magnitude).
                normal_stiffness: Normal stiffness (used for normal force magnitude).
                dmin: Minimum distance (used for normal force magnitude).
                no_mu: Whether to not multiply by mu

            Returns:
                Friction force
            )ipc_Qu8mg5v7",
            "collision"_a, "rest_positions"_a, "lagged_displacements"_a,
            "velocities"_a, "normal_potential"_a, "normal_stiffness"_a,
            "dmin"_a = 0, "no_mu"_a = false)
        .def(
            "force_jacobian",
            py::overload_cast<
                const TangentialCollision&, Eigen::ConstRef<VectorMax12d>,
                Eigen::ConstRef<VectorMax12d>, Eigen::ConstRef<VectorMax12d>,
                const NormalPotential&, const double,
                const TangentialPotential::DiffWRT, const double>(
                &TangentialPotential::force_jacobian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the friction force Jacobian.

            Parameters:
                collision: The collision
                rest_positions: Rest positions of the vertices (rowwise).
                lagged_displacements: Previous displacements of the vertices (rowwise).
                velocities: Current displacements of the vertices (rowwise).
                normal_potential: Normal potential (used for normal force magnitude).
                normal_stiffness: Normal stiffness (used for normal force magnitude).
                wrt: Variable to differentiate the friction force with respect to.
                dmin: Minimum distance (used for normal force magnitude).

            Returns:
                Friction force Jacobian
            )ipc_Qu8mg5v7",
            "collision"_a, "rest_positions"_a, "lagged_displacements"_a,
            "velocities"_a, "normal_potential"_a, "normal_stiffness"_a, "wrt"_a,
            "dmin"_a = 0);
}
