#include <common.hpp>

#include <ipc/friction/friction_collisions.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_collisions(py::module_& m)
{
    py::class_<FrictionCollisions>(m, "FrictionCollisions")
        .def(py::init())
        // Original build method for a single friction coefficient
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const Collisions&,
                const BarrierPotential&, double, double>(
                &FrictionCollisions::build),
            py::arg("mesh"), py::arg("vertices"), py::arg("collisions"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("mu"))
        // Overload for build method with a vector of friction coefficients
        .def(
            "build",
            [](FrictionCollisions& self, const CollisionMesh& mesh,
               const Eigen::MatrixXd& vertices, const Collisions& collisions,
               const BarrierPotential& barrier_potential,
               const double barrier_stiffness, const Eigen::VectorXd& mus) {
                self.build(
                    mesh, vertices, collisions, barrier_potential,
                    barrier_stiffness, mus);
            },
            py::arg("mesh"), py::arg("vertices"), py::arg("collisions"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("mus"))
        // Overload for build method with a blend function for friction coefficients
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const Collisions&,
                const BarrierPotential&, const double, const Eigen::VectorXd&,
                const std::function<double(double, double)>&>(
                &FrictionCollisions::build),
            py::arg("mesh"), py::arg("vertices"), py::arg("collisions"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("mus"), py::arg("blend_mu"))
        // New overload for build method supporting static/kinetic and pairwise friction coefficients
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const Collisions&,
                const BarrierPotential&, const double, const double, const double,
                const std::map<std::tuple<int, int>, std::pair<double, double>>&>(
                &FrictionCollisions::build),
            py::arg("mesh"), py::arg("vertices"), py::arg("collisions"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("static_mu"), py::arg("kinetic_mu"),
            py::arg("pairwise_friction"),
            R"ipc_Qu8mg5v7(
            Build the friction collisions with static, kinetic, and pairwise friction coefficients.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertex positions.
                collisions: Collision data.
                barrier_potential: The barrier potential used for normal force.
                barrier_stiffness: Barrier stiffness.
                static_mu: Global static friction coefficient.
                kinetic_mu: Global kinetic friction coefficient.
                pairwise_friction: Pairwise static and kinetic friction coefficients.
            )ipc_Qu8mg5v7")
        // Function to get the number of friction collisions
        .def("__len__", &FrictionCollisions::size,
            "Get the number of friction collisions.")
        // Check if the friction collisions are empty
        .def("empty", &FrictionCollisions::empty,
            "Get if the friction collisions are empty.")
        // Clear the friction collisions
        .def("clear", &FrictionCollisions::clear,
            "Clear the friction collisions.")
        // Access an individual collision at index i
        .def(
            "__getitem__",
            [](FrictionCollisions& self, size_t i) -> FrictionCollision& {
                return self[i];
            },
            py::return_value_policy::reference,
            R"ipc_Qu8mg5v7(
            Get a reference to collision at index i.

            Parameters:
                i: The index of the collision.

            Returns:
                A reference to the collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        // Default blend function for combining friction coefficients
        .def_static("default_blend_mu", &FrictionCollisions::default_blend_mu,
            py::arg("mu0"), py::arg("mu1"))
        // Expose the different types of friction collisions for external manipulation
        .def_readwrite("vv_collisions", &FrictionCollisions::vv_collisions)
        .def_readwrite("ev_collisions", &FrictionCollisions::ev_collisions)
        .def_readwrite("ee_collisions", &FrictionCollisions::ee_collisions)
        .def_readwrite("fv_collisions", &FrictionCollisions::fv_collisions);
}
