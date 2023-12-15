#include <common.hpp>

#include <ipc/friction/friction_collisions.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_collisions(py::module_& m)
{
    py::class_<FrictionCollisions>(m, "FrictionCollisions")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const Collisions&,
                const BarrierPotential&, double, double>(
                &FrictionCollisions::build),
            py::arg("mesh"), py::arg("vertices"), py::arg("collisions"),
            py::arg("barrier_potential"), py::arg("barrier_stiffness"),
            py::arg("mu"))
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
        .def(
            "__len__", &FrictionCollisions::size,
            "Get the number of friction collisions.")
        .def(
            "empty", &FrictionCollisions::empty,
            "Get if the friction collisions are empty.")
        .def(
            "clear", &FrictionCollisions::clear,
            "Clear the friction collisions.")
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
        .def_static(
            "default_blend_mu", &FrictionCollisions::default_blend_mu,
            py::arg("mu0"), py::arg("mu1"))
        .def_readwrite("vv_collisions", &FrictionCollisions::vv_collisions)
        .def_readwrite("ev_collisions", &FrictionCollisions::ev_collisions)
        .def_readwrite("ee_collisions", &FrictionCollisions::ee_collisions)
        .def_readwrite("fv_collisions", &FrictionCollisions::fv_collisions);
}
