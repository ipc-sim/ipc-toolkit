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
                const NormalCollisions&, const NormalPotential&, double>(
                &TangentialCollisions::build),
            "mesh"_a, "vertices"_a, "collisions"_a, "normal_potential"_a,
            "mu"_a)
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                const NormalCollisions&, const NormalPotential&, double,
                double>(&TangentialCollisions::build),
            "mesh"_a, "vertices"_a, "collisions"_a, "normal_potential"_a,
            "mu_s"_a, "mu_k"_a)
        .def(
            "build",
            [](TangentialCollisions& self, const CollisionMesh& mesh,
               Eigen::ConstRef<Eigen::MatrixXd> vertices,
               const NormalCollisions& collisions,
               const NormalPotential& normal_potential,
               Eigen::ConstRef<Eigen::VectorXd> mu_s,
               Eigen::ConstRef<Eigen::VectorXd> mu_k) {
                self.build(
                    mesh, vertices, collisions, normal_potential, mu_s, mu_k);
            },
            "mesh"_a, "vertices"_a, "collisions"_a, "normal_potential"_a,
            "mu_s"_a, "mu_k"_a)
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                const NormalCollisions&, const NormalPotential&,
                Eigen::ConstRef<Eigen::VectorXd>,
                Eigen::ConstRef<Eigen::VectorXd>,
                const std::function<double(double, double)>&>(
                &TangentialCollisions::build),
            "mesh"_a, "vertices"_a, "collisions"_a, "normal_potential"_a,
            "mu_s"_a, "mu_k"_a, "blend_mu"_a)
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
            "reset_lagged_anisotropic_friction_coefficients",
            &TangentialCollisions::reset_lagged_anisotropic_friction_coefficients,
            "Set lagged effective μ to scalar mu_s/mu_k on each collision (done automatically after build).")
        .def(
            "update_lagged_anisotropic_friction_coefficients",
            &TangentialCollisions::update_lagged_anisotropic_friction_coefficients,
            "mesh"_a, "rest_positions"_a, "lagged_displacements"_a, "velocities"_a,
            "Refresh matchstick effective μ from lagged geometry and slip. "
            "Call when mu_s_aniso is nonzero (e.g. each Newton iteration).")
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
