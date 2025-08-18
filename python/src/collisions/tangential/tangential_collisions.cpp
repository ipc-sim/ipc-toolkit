#include <common.hpp>

#include <ipc/collisions/tangential/tangential_collisions.hpp>

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
            "mesh"_a, "vertices"_a, "collisions"_a, "normal_potential"_a,
            "normal_stiffness"_a, "mu"_a)
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                const NormalCollisions&, const NormalPotential&, double, double,
                double>(&TangentialCollisions::build),
            "mesh"_a, "vertices"_a, "collisions"_a, "normal_potential"_a,
            "normal_stiffness"_a, "mu_s"_a, "mu_k"_a)
        .def(
            "build",
            [](TangentialCollisions& self, const CollisionMesh& mesh,
               Eigen::ConstRef<Eigen::MatrixXd> vertices,
               const NormalCollisions& collisions,
               const NormalPotential& normal_potential,
               const double normal_stiffness,
               Eigen::ConstRef<Eigen::VectorXd> mu_s,
               Eigen::ConstRef<Eigen::VectorXd> mu_k) {
                self.build(
                    mesh, vertices, collisions, normal_potential,
                    normal_stiffness, mu_s, mu_k);
            },
            "mesh"_a, "vertices"_a, "collisions"_a, "normal_potential"_a,
            "normal_stiffness"_a, "mu_s"_a, "mu_k"_a)
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                const NormalCollisions&, const NormalPotential&, const double,
                Eigen::ConstRef<Eigen::VectorXd>,
                Eigen::ConstRef<Eigen::VectorXd>,
                const std::function<double(double, double)>&>(
                &TangentialCollisions::build),
            "mesh"_a, "vertices"_a, "collisions"_a, "normal_potential"_a,
            "normal_stiffness"_a, "mu_s"_a, "mu_k"_a, "blend_mu"_a)
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
            "i"_a)
        .def_static(
            "default_blend_mu", &TangentialCollisions::default_blend_mu,
            "mu0"_a, "mu1"_a)
        .def_readwrite("vv_collisions", &TangentialCollisions::vv_collisions)
        .def_readwrite("ev_collisions", &TangentialCollisions::ev_collisions)
        .def_readwrite("ee_collisions", &TangentialCollisions::ee_collisions)
        .def_readwrite("fv_collisions", &TangentialCollisions::fv_collisions);
}
