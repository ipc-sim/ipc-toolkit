#include <common.hpp>

#include <ipc/collisions/normal/normal_collisions.hpp>

namespace py = pybind11;
using namespace ipc;

void define_normal_collisions(py::module_& m)
{
    py::class_<NormalCollisions>(m, "NormalCollisions")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const double,
                const double, const BroadPhaseMethod>(&NormalCollisions::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of collisions used to compute the barrier potential.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.
                dmin: Minimum distance.
                broad_phase_method: Broad-phase method to use.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("dhat"),
            py::arg("dmin") = 0,
            py::arg("broad_phase_method") = DEFAULT_BROAD_PHASE_METHOD)
        .def(
            "build",
            py::overload_cast<
                const Candidates&, const CollisionMesh&, const Eigen::MatrixXd&,
                const double, const double>(&NormalCollisions::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of collisions used to compute the barrier potential.

            Parameters:
                candidates: Distance candidates from which the collision set is built.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.
                dmin:  Minimum distance.
            )ipc_Qu8mg5v7",
            py::arg("candidates"), py::arg("mesh"), py::arg("vertices"),
            py::arg("dhat"), py::arg("dmin") = 0)
        .def(
            "compute_minimum_distance",
            &NormalCollisions::compute_minimum_distance,
            R"ipc_Qu8mg5v7(
            Computes the minimum distance between any non-adjacent elements.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The minimum distance between any non-adjacent elements.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"))
        .def(
            "__len__", &NormalCollisions::size, "Get the number of collisions.")
        .def(
            "empty", &NormalCollisions::empty,
            "Get if the collision set are empty.")
        .def("clear", &NormalCollisions::clear, "Clear the collision set.")
        .def(
            "__getitem__",
            [](NormalCollisions& self, size_t i) -> NormalCollision& {
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
        .def(
            "is_vertex_vertex", &NormalCollisions::is_vertex_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is a vertex-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is a vertex-vertex collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "is_edge_vertex", &NormalCollisions::is_edge_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an edge-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an edge-vertex collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "is_edge_edge", &NormalCollisions::is_edge_edge,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an edge-edge collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an edge-edge collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "is_face_vertex", &NormalCollisions::is_face_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an face-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an face-vertex collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "is_plane_vertex", &NormalCollisions::is_plane_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an plane-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an plane-vertex collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "to_string", &NormalCollisions::to_string, py::arg("mesh"),
            py::arg("vertices"))
        .def_property(
            "use_convergent_formulation",
            &NormalCollisions::use_convergent_formulation,
            &NormalCollisions::set_use_convergent_formulation,
            "If the collisions should use the convergent formulation.")
        .def_property(
            "are_shape_derivatives_enabled",
            &NormalCollisions::are_shape_derivatives_enabled,
            &NormalCollisions::set_are_shape_derivatives_enabled,
            "If the collisions are using the convergent formulation.")
        .def_readwrite("vv_collisions", &NormalCollisions::vv_collisions)
        .def_readwrite("ev_collisions", &NormalCollisions::ev_collisions)
        .def_readwrite("ee_collisions", &NormalCollisions::ee_collisions)
        .def_readwrite("fv_collisions", &NormalCollisions::fv_collisions)
        .def_readwrite("pv_collisions", &NormalCollisions::pv_collisions);
}
