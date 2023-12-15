#include <common.hpp>

#include <ipc/collisions/collisions.hpp>

namespace py = pybind11;
using namespace ipc;

void define_collisions(py::module_& m)
{
    py::class_<Collisions>(m, "Collisions")
        .def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const double,
                const double, const BroadPhaseMethod>(&Collisions::build),
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
                const double, const double>(&Collisions::build),
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
            "compute_minimum_distance", &Collisions::compute_minimum_distance,
            R"ipc_Qu8mg5v7(
            Computes the minimum distance between any non-adjacent elements.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The minimum distance between any non-adjacent elements.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"))
        .def("__len__", &Collisions::size, "Get the number of collisions.")
        .def("empty", &Collisions::empty, "Get if the collision set is empty.")
        .def("clear", &Collisions::clear, "Clear the collision set.")
        .def(
            "__getitem__",
            [](Collisions& self, size_t i) -> Collision& { return self[i]; },
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
            "is_vertex_vertex", &Collisions::is_vertex_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is a vertex-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is a vertex-vertex collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "is_edge_vertex", &Collisions::is_edge_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an edge-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an edge-vertex collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "is_edge_edge", &Collisions::is_edge_edge,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an edge-edge collision.

            Parameters:
                i: The index of the collision.nose

            Returns:
                If the collision at i is an edge-edge collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "is_face_vertex", &Collisions::is_face_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an face-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an face-vertex collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "is_plane_vertex", &Collisions::is_plane_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an plane-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an plane-vertex collision.
            )ipc_Qu8mg5v7",
            py::arg("i"))
        .def(
            "to_string", &Collisions::to_string, "", py::arg("mesh"),
            py::arg("vertices"))
        .def_property(
            "use_convergent_formulation",
            &Collisions::use_convergent_formulation,
            &Collisions::set_use_convergent_formulation,
            "If the collisions should use the convergent formulation.")
        .def_property(
            "are_shape_derivatives_enabled",
            &Collisions::are_shape_derivatives_enabled,
            &Collisions::set_are_shape_derivatives_enabled,
            "If the collisions are using the convergent formulation.")
        .def(
            "to_string", &Collisions::to_string, py::arg("mesh"),
            py::arg("vertices"))
        .def_readwrite("vv_collisions", &Collisions::vv_collisions)
        .def_readwrite("ev_collisions", &Collisions::ev_collisions)
        .def_readwrite("ee_collisions", &Collisions::ee_collisions)
        .def_readwrite("fv_collisions", &Collisions::fv_collisions)
        .def_readwrite("pv_collisions", &Collisions::pv_collisions);
}
