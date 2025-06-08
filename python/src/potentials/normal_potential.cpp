#include <common.hpp>
#include <potentials/potential.hpp>

#include <ipc/potentials/normal_potential.hpp>

namespace py = pybind11;
using namespace ipc;

void define_normal_potential(py::module_& m)
{
    py::class_<NormalPotential> normal_potential(m, "NormalPotential");

    define_potential_methods<NormalCollisions>(normal_potential);

    normal_potential
        .def(
            "shape_derivative",
            py::overload_cast<
                const NormalCollisions&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>>(
                &NormalPotential::shape_derivative, py::const_),
            R"ipc_Qu8mg5v7(
            Compute the shape derivative of the potential.

            std::runtime_error If the collision collisions were not built with shape derivatives enabled.

            Parameters:
                collisions: The set of collisions.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The derivative of the force with respect to X, the rest vertices.
            )ipc_Qu8mg5v7",
            py::arg("collisions"), py::arg("mesh"), py::arg("vertices"))
        .def(
            "shape_derivative",
            [](const NormalPotential& self, const NormalCollision& collision,
               const std::array<index_t, 4>& vertex_ids,
               Eigen::ConstRef<VectorMax12d> rest_positions,
               Eigen::ConstRef<VectorMax12d> positions) {
                std::vector<Eigen::Triplet<double>> out;
                self.shape_derivative(
                    collision, vertex_ids, rest_positions, positions, out);
                return out;
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
        .def(
            "force_magnitude", &NormalPotential::force_magnitude,
            R"ipc_Qu8mg5v7(
            Compute the force magnitude for a collision.

            Parameters:
                distance_squared: The squared distance between elements.
                dmin: The minimum distance offset to the barrier.
                barrier_stiffness: The barrier stiffness.

            Returns:
                The force magnitude.
            )ipc_Qu8mg5v7",
            py::arg("distance_squared"), py::arg("dmin"),
            py::arg("barrier_stiffness"))
        .def(
            "force_magnitude_gradient",
            &NormalPotential::force_magnitude_gradient,
            R"ipc_Qu8mg5v7(
            Compute the gradient of the force magnitude for a collision.

            Parameters:
                distance_squared: The squared distance between elements.
                distance_squared_gradient: The gradient of the squared distance.
                dmin: The minimum distance offset to the barrier.
                barrier_stiffness: The stiffness of the barrier.

            Returns:
                The gradient of the force.
            )ipc_Qu8mg5v7",
            py::arg("distance_squared"), py::arg("distance_squared_gradient"),
            py::arg("dmin"), py::arg("barrier_stiffness"));
}
