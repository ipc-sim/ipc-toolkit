#include <common.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/smooth_contact/smooth_contact_potential.hpp>

using namespace ipc;

void define_barrier_potential(py::module_& m)
{
    py::class_<BarrierPotential, NormalPotential>(m, "BarrierPotential")
        .def(
            py::init<const double, const bool>(),
            R"ipc_Qu8mg5v7(
            Construct a barrier potential.

            Parameters:
                dhat: The activation distance of the barrier.
            )ipc_Qu8mg5v7",
            "dhat"_a, "use_physical_barrier"_a = false)
        .def(
            py::init<
                const std::shared_ptr<Barrier>, const double, const bool>(),
            R"ipc_Qu8mg5v7(
            Construct a barrier potential.

            Parameters:
                barrier: The barrier function.
                dhat: The activation distance of the barrier.
            )ipc_Qu8mg5v7",
            "barrier"_a, "dhat"_a, "use_physical_barrier"_a = false)
        .def_property(
            "dhat", &BarrierPotential::dhat, &BarrierPotential::set_dhat,
            "Barrier activation distance.")
        .def_property(
            "barrier",
            py::cpp_function(
                &BarrierPotential::barrier, py::return_value_policy::reference),
            &BarrierPotential::set_barrier,
            "Barrier function used to compute the potential.");
}

void define_smooth_potential(py::module_& m)
{
    py::class_<SmoothContactParameters>(m, "SmoothContactParameters")
        .def(
            py::init<
                const double, const double, const double, const double,
                const double, const int>(),
            R"ipc_Qu8mg5v7(
            Construct parameter set for smooth contact.

            Parameters:
                dhat, alpha_t, beta_t, alpha_n, beta_n, r
            )ipc_Qu8mg5v7",
            "dhat"_a, "alpha_t"_a, "beta_t"_a, "alpha_n"_a, "beta_n"_a, "r"_a)
        .def(
            py::init<const double, const double, const double, const int>(),
            R"ipc_Qu8mg5v7(
            Construct parameter set for smooth contact.

            Parameters:
                dhat, alpha_t, beta_t, r
            )ipc_Qu8mg5v7",
            "dhat"_a, "alpha_t"_a, "beta_t"_a, "r"_a)
        .def_readonly("dhat", &SmoothContactParameters::dhat)
        .def_readonly("alpha_t", &SmoothContactParameters::alpha_t)
        .def_readonly("beta_t", &SmoothContactParameters::beta_t)
        .def_readonly("alpha_n", &SmoothContactParameters::alpha_n)
        .def_readonly("beta_n", &SmoothContactParameters::beta_n)
        .def_readonly("r", &SmoothContactParameters::r);

    py::class_<SmoothContactPotential>(m, "SmoothPotential")
        .def(
            py::init<const SmoothContactParameters&>(),
            R"ipc_Qu8mg5v7(
            Construct a smooth barrier potential.

            Parameters:
                param: A set of parameters.
            )ipc_Qu8mg5v7",
            "param"_a)
        .def(
            "__call__",
            py::overload_cast<
                const SmoothCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &ipc::SmoothContactPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the barrier potential for a set of collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The sum of all barrier potentials (not scaled by the barrier stiffness).
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a)
        .def(
            "gradient",
            py::overload_cast<
                const SmoothCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &ipc::SmoothContactPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a)
        .def(
            "hessian",
            py::overload_cast<
                const SmoothCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>, const PSDProjectionMethod>(
                &ipc::SmoothContactPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                project_hessian_to_psd: Make sure the hessian is positive semi-definite.

            Returns:
                The hessian of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|x|vertices|.
            )ipc_Qu8mg5v7",
            "collisions"_a, "mesh"_a, "vertices"_a,
            "project_hessian_to_psd"_a = PSDProjectionMethod::NONE)
        .def(
            "__call__",
            py::overload_cast<
                const SmoothCollision&, Eigen::ConstRef<Eigen::VectorXd>>(
                &ipc::SmoothContactPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a)
        .def(
            "gradient",
            py::overload_cast<
                const SmoothCollision&, Eigen::ConstRef<Eigen::VectorXd>>(
                &SmoothContactPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The gradient of the potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a)
        .def(
            "hessian",
            py::overload_cast<
                const SmoothCollision&, Eigen::ConstRef<Eigen::VectorXd>,
                const PSDProjectionMethod>(
                &SmoothContactPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The hessian of the potential.
            )ipc_Qu8mg5v7",
            "collision"_a, "x"_a,
            "project_hessian_to_psd"_a = PSDProjectionMethod::NONE);
}
