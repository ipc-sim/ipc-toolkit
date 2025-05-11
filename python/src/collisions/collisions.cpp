#include <common.hpp>

#include <ipc/collisions/collisions.hpp>
#include <ipc/smooth_contact/smooth_collisions.hpp>

namespace py = pybind11;
using namespace ipc;

template <typename collision_type, typename parent_type>
void define_SmoothCollisionTemplate(py::module_& m, std::string name)
{
    py::class_<collision_type, parent_type>(m, name.c_str())
    .def("name", &collision_type::name, "Get the type name of collision")
    .def("num_vertices", &collision_type::num_vertices, "Get the number of vertices");
}

template <int dim>
void define_SmoothCollisions(py::module_& m, std::string name)
{
    py::class_<SmoothCollisions<dim>>(m, name.c_str())
        .def(py::init())
        .def_readwrite("use_high_order_quadrature", &SmoothCollisions<dim>::use_high_order_quadrature)
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, const Eigen::MatrixXd&, const ParameterType, const bool, const BroadPhaseMethod>(&SmoothCollisions<dim>::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of collisions used to compute the barrier potential.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                param: ParameterType.
                use_adaptive_dhat: If the adaptive dhat should be used.
                broad_phase: Broad phase method.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"), py::arg("param"),
            py::arg("use_adaptive_dhat") = false, py::arg("broad_phase") = DEFAULT_BROAD_PHASE_METHOD)
        .def(
            "compute_minimum_distance", &SmoothCollisions<dim>::compute_minimum_distance,
            R"ipc_Qu8mg5v7(
            Computes the minimum distance between any non-adjacent elements.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.

            Returns:
                The minimum distance between any non-adjacent elements.
            )ipc_Qu8mg5v7",
            py::arg("mesh"), py::arg("vertices"))
        .def("__len__", &SmoothCollisions<dim>::size, "Get the number of collisions.")
        .def("empty", &SmoothCollisions<dim>::empty, "Get if the collision set is empty.")
        .def("clear", &SmoothCollisions<dim>::clear, "Clear the collision set.")
        .def(
            "__getitem__",
            [](SmoothCollisions<dim>& self, size_t i) -> typename SmoothCollisions<dim>::value_type& { return self[i]; },
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
            "to_string", &SmoothCollisions<dim>::to_string, py::arg("mesh"),
            py::arg("vertices"), py::arg("param"))
        .def("n_candidates", &SmoothCollisions<dim>::n_candidates, "Get the number of candidates.");
}

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
            [](Collisions& self, size_t i) -> Collision<4>& { return self[i]; },
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

    py::class_<SmoothCollision<6>>(m, "SmoothCollision2")
        .def("n_dofs", &SmoothCollision<6>::n_dofs, "Get the degree of freedom")
        .def(
            "__call__",
            &SmoothCollision<6>::operator(),
            R"ipc_Qu8mg5v7(
            Compute the potential.
    
            Parameters:
                positions: The vertex positions.
                params: The parameters.
    
            Returns:
                The potential (not scaled by the barrier stiffness) of this collision pair.
            )ipc_Qu8mg5v7",
            py::arg("positions"), py::arg("params"))
        .def(
            "__getitem__",
            [](SmoothCollision<6>& self, size_t i) -> long { return self[i]; },
            R"ipc_Qu8mg5v7(
            Get primitive id.

            Parameters:
                i: 0 or 1.

            Returns:
                The index of the primitive.
            )ipc_Qu8mg5v7",
            py::arg("i"));

    define_SmoothCollisionTemplate<SmoothCollisionTemplate<max_vert_2d, Edge2, Point2>, SmoothCollision<6>>(m, "Edge2Point2Collision");
    define_SmoothCollisionTemplate<SmoothCollisionTemplate<max_vert_2d, Point2, Point2>, SmoothCollision<6>>(m, "Point2Point2Collision");

    define_SmoothCollisions<2>(m, "SmoothCollisions2");
    define_SmoothCollisions<3>(m, "SmoothCollisions3");
}
