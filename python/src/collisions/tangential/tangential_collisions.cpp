#include <common.hpp>

#include <ipc/collisions/tangential/tangential_collisions.hpp>

namespace py = pybind11;
using namespace ipc;

void define_tangential_collisions(py::module_& m)
{
    py::class_<TangentialCollisions>(m, "TangentialCollisions")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                const NormalCollisions&, const NormalPotential&, double,
                double>(&TangentialCollisions::build),
            py::arg("mesh"), py::arg("vertices"), py::arg("collisions"),
            py::arg("normal_potential"), py::arg("normal_stiffness"),
            py::arg("mu"))
        .def(
            "build",
            [](TangentialCollisions& self, const CollisionMesh& mesh,
               Eigen::ConstRef<Eigen::MatrixXd> vertices,
               const NormalCollisions& collisions,
               const NormalPotential& normal_potential,
               const double normal_stiffness,
               Eigen::ConstRef<Eigen::VectorXd> mus) {
                self.build(
                    mesh, vertices, collisions, normal_potential,
                    normal_stiffness, mus);
            },
            py::arg("mesh"), py::arg("vertices"), py::arg("collisions"),
            py::arg("normal_potential"), py::arg("normal_stiffness"),
            py::arg("mus"))
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                const NormalCollisions&, const NormalPotential&, const double,
                Eigen::ConstRef<Eigen::VectorXd>,
                const std::function<double(double, double)>&>(
                &TangentialCollisions::build),
            py::arg("mesh"), py::arg("vertices"), py::arg("collisions"),
            py::arg("normal_potential"), py::arg("normal_stiffness"),
            py::arg("mus"), py::arg("blend_mu"))
        .def(
            "__len__", &TangentialCollisions::size,
            "Get the number of friction collisions.")
        .def(
            "empty", &TangentialCollisions::empty,
            "Get if the friction collisions are empty.")
        .def(
            "clear", &TangentialCollisions::clear,
            "Clear the friction collisions.")
        .def(
            "__getitem__",
            [](TangentialCollisions& self, size_t i) -> TangentialCollision& {
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
            "default_blend_mu", &TangentialCollisions::default_blend_mu,
            py::arg("mu0"), py::arg("mu1"))
        .def_readwrite("vv_collisions", &TangentialCollisions::vv_collisions)
        .def_readwrite("ev_collisions", &TangentialCollisions::ev_collisions)
        .def_readwrite("ee_collisions", &TangentialCollisions::ee_collisions)
        .def_readwrite("fv_collisions", &TangentialCollisions::fv_collisions);
}
