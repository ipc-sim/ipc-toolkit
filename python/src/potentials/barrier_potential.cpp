#include <common.hpp>

#include <ipc/potentials/barrier_potential.hpp>
#include <ipc/smooth_contact/smooth_contact_potential.hpp>

namespace py = pybind11;
using namespace ipc;

void define_barrier_potential(py::module_& m)
{
    py::class_<BarrierPotential>(m, "BarrierPotential")
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a barrier potential.

            Parameters:
                dhat: The activation distance of the barrier.
            )ipc_Qu8mg5v7",
            py::arg("dhat"))
        .def(
            "__call__",
            py::overload_cast<
                const Collisions&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
                &BarrierPotential::Potential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the barrier potential for a set of collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The sum of all barrier potentials (not scaled by the barrier stiffness).
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "gradient",
            py::overload_cast<
                const Collisions&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
                &BarrierPotential::Potential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "hessian",
            py::overload_cast<
                const Collisions&, const CollisionMesh&, const Eigen::MatrixXd&,
                const bool>(&BarrierPotential::Potential::hessian, py::const_),
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
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "shape_derivative",
            py::overload_cast<
                const Collisions&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
                &BarrierPotential::shape_derivative, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the shape derivative of the potential.

            std::runtime_error If the collision collisions were not built
            with shape derivatives enabled.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The derivative of the force with respect to X, the rest
                vertices.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "__call__",
            py::overload_cast<const Collision<4>&, const VectorMax12d&>(
                &BarrierPotential::operator(), py::const_),
            R"ipc_Qu8mg5v7(
            Compute the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"))
        .def(
            "gradient",
            py::overload_cast<const Collision<4>&, const VectorMax12d&>(
                &BarrierPotential::gradient, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The gradient of the potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"))
        .def(
            "hessian",
            py::overload_cast<
                const Collision<4>&, const VectorMax12d&, const bool>(
                &BarrierPotential::hessian, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The hessian of the potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "shape_derivative",
            [](const DistanceBasedPotential& self,
               const Collision<4>& collision,
               const std::array<long, 4>& vertex_ids,
               const VectorMax12d& rest_positions,
               const VectorMax12d& positions) {
                std::vector<Eigen::Triplet<double>> out;
                self.shape_derivative(
                    collision, vertex_ids, rest_positions, positions, out);
                Eigen::SparseMatrix<double> out_mat(
                    rest_positions.size(), rest_positions.size());
                out_mat.setFromTriplets(out.begin(), out.end());
                return out_mat;
            },
            R"ipc_Qu8mg5v7(
            Compute the shape derivative of the potential for a single collision.

            Parameters:
                collision: The collision.
                vertex_ids: The collision stencil's vertex ids.
                rest_positions: The collision stencil's rest positions.
                positions: The collision stencil's positions.
                ,out]: out Store the triplets of the shape derivative here.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("vertex_ids"),
            py::arg("rest_positions"), py::arg("positions"))
        .def_property(
            "dhat", &BarrierPotential::dhat, &BarrierPotential::set_dhat,
            "Barrier activation distance.")
        .def_property(
            "barrier",
            py::cpp_function(
                &BarrierPotential::barrier, py::return_value_policy::reference),
            &BarrierPotential::set_barrier, "Barrier activation distance.");
}

void define_smooth_potential(py::module_& m)
{
    py::class_<ParameterType>(m, "ParameterType")
        .def(
            py::init<
                const double&, const double&, const double&, const double&,
                const double&, const int&>(),
            R"ipc_Qu8mg5v7(
            Construct parameter set for smooth contact.

            Parameters:
                dhat, alpha_t, beta_t, alpha_n, beta_n, r
            )ipc_Qu8mg5v7",
            py::arg("dhat"), py::arg("alpha_t"), py::arg("beta_t"),
            py::arg("alpha_n"), py::arg("beta_n"), py::arg("r"))
        .def(
            py::init<const double&, const double&, const double&, const int&>(),
            R"ipc_Qu8mg5v7(
            Construct parameter set for smooth contact.

            Parameters:
                dhat, alpha_t, beta_t, r
            )ipc_Qu8mg5v7",
            py::arg("dhat"), py::arg("alpha_t"), py::arg("beta_t"),
            py::arg("r"))
        .def_readonly("dhat", &ParameterType::dhat)
        .def_readonly("alpha_t", &ParameterType::alpha_t)
        .def_readonly("beta_t", &ParameterType::beta_t)
        .def_readonly("alpha_n", &ParameterType::alpha_n)
        .def_readonly("beta_n", &ParameterType::beta_n)
        .def_readonly("r", &ParameterType::r);

    py::class_<SmoothContactPotential<SmoothCollisions<2>>>(
        m, "SmoothPotential")
        .def(
            py::init<const ParameterType&>(),
            R"ipc_Qu8mg5v7(
            Construct a smooth barrier potential.

            Parameters:
                param: A set of parameters.
            )ipc_Qu8mg5v7",
            py::arg("param"))
        .def(
            "__call__",
            py::overload_cast<
                const SmoothCollisions<2>&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
                &SmoothContactPotential<
                    SmoothCollisions<2>>::Potential::operator(),
                py::const_),
            R"ipc_Qu8mg5v7(
            Compute the barrier potential for a set of collisions.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The sum of all barrier potentials (not scaled by the barrier stiffness).
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "gradient",
            py::overload_cast<
                const SmoothCollisions<2>&, const CollisionMesh&,
                const Eigen::MatrixXd&>(
                &SmoothContactPotential<
                    SmoothCollisions<2>>::Potential::gradient,
                py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the barrier potential.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The gradient of all barrier potentials (not scaled by the barrier stiffness). This will have a size of |vertices|.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "hessian",
            py::overload_cast<
                const SmoothCollisions<2>&, const CollisionMesh&,
                const Eigen::MatrixXd&, const bool>(
                &SmoothContactPotential<
                    SmoothCollisions<2>>::Potential::hessian,
                py::const_),
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
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"),
            py::arg("project_hessian_to_psd") = false)
        .def(
            "__call__",
            py::overload_cast<
                const SmoothCollision<6>&, const Vector<double, -1, 18>&>(
                &SmoothContactPotential<SmoothCollisions<2>>::operator(),
                py::const_),
            R"ipc_Qu8mg5v7(
            Compute the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"))
        .def(
            "gradient",
            py::overload_cast<
                const SmoothCollision<6>&, const Vector<double, -1, 18>&>(
                &SmoothContactPotential<SmoothCollisions<2>>::gradient,
                py::const_),
            R"ipc_Qu8mg5v7(
            Compute the gradient of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The gradient of the potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"))
        .def(
            "hessian",
            py::overload_cast<
                const SmoothCollision<6>&, const Vector<double, -1, 18>&,
                const bool>(
                &SmoothContactPotential<SmoothCollisions<2>>::hessian,
                py::const_),
            R"ipc_Qu8mg5v7(
            Compute the hessian of the potential for a single collision.

            Parameters:
                collision: The collision.
                x: The collision stencil's degrees of freedom.

            Returns:
                The hessian of the potential.
            )ipc_Qu8mg5v7",
            py::arg("collision"), py::arg("x"),
            py::arg("project_hessian_to_psd") = false);
}
