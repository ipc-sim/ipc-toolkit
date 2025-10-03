#include <common.hpp>

#include <ipc/collisions/normal/normal_collisions.hpp>

using namespace ipc;

void define_normal_collisions(py::module_& m)
{
    py::class_<NormalCollisions> normal_collisions(m, "NormalCollisions");

    py::enum_<NormalCollisions::CollisionSetType>(
        normal_collisions, "CollisionSetType")
        .value(
            "IPC", NormalCollisions::CollisionSetType::IPC,
            "IPC set type, which uses the original formulation described in [Li et al. 2020].")
        .value(
            "IMPROVED_MAX_APPROX",
            NormalCollisions::CollisionSetType::IMPROVED_MAX_APPROX,
            "Improved max approximation set type, which uses the improved max approximation formulation described in [Li et al. 2023].")
        .value(
            "OGC", NormalCollisions::CollisionSetType::OGC,
            "Offset Geometric Contact set type, which uses the formulation described in [Chen et al. 2025].")
        .export_values();

    normal_collisions.def(py::init())
        .def(
            "build",
            py::overload_cast<
                const CollisionMesh&, Eigen::ConstRef<Eigen::MatrixXd>,
                const double, const double, const std::shared_ptr<BroadPhase>&>(
                &NormalCollisions::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of collisions used to compute the barrier potential.

            Parameters:
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.
                dmin: Minimum distance.
                broad_phase: Broad-phase to use.
            )ipc_Qu8mg5v7",
            "mesh"_a, "vertices"_a, "dhat"_a, "dmin"_a = 0,
            "broad_phase"_a = make_default_broad_phase())
        .def(
            "build",
            py::overload_cast<
                const Candidates&, const CollisionMesh&,
                Eigen::ConstRef<Eigen::MatrixXd>, const double, const double>(
                &NormalCollisions::build),
            R"ipc_Qu8mg5v7(
            Initialize the set of collisions used to compute the barrier potential.

            Parameters:
                candidates: Distance candidates from which the collision set is built.
                mesh: The collision mesh.
                vertices: Vertices of the collision mesh.
                dhat: The activation distance of the barrier.
                dmin:  Minimum distance.
            )ipc_Qu8mg5v7",
            "candidates"_a, "mesh"_a, "vertices"_a, "dhat"_a, "dmin"_a = 0)
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
            "mesh"_a, "vertices"_a)
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
            "i"_a)
        .def(
            "is_vertex_vertex", &NormalCollisions::is_vertex_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is a vertex-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is a vertex-vertex collision.
            )ipc_Qu8mg5v7",
            "i"_a)
        .def(
            "is_edge_vertex", &NormalCollisions::is_edge_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an edge-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an edge-vertex collision.
            )ipc_Qu8mg5v7",
            "i"_a)
        .def(
            "is_edge_edge", &NormalCollisions::is_edge_edge,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an edge-edge collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an edge-edge collision.
            )ipc_Qu8mg5v7",
            "i"_a)
        .def(
            "is_face_vertex", &NormalCollisions::is_face_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an face-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an face-vertex collision.
            )ipc_Qu8mg5v7",
            "i"_a)
        .def(
            "is_plane_vertex", &NormalCollisions::is_plane_vertex,
            R"ipc_Qu8mg5v7(
            Get if the collision at i is an plane-vertex collision.

            Parameters:
                i: The index of the collision.

            Returns:
                If the collision at i is an plane-vertex collision.
            )ipc_Qu8mg5v7",
            "i"_a)
        .def("__str__", &NormalCollisions::to_string, "mesh"_a, "vertices"_a)
        .def_property(
            "use_area_weighting", &NormalCollisions::use_area_weighting,
            &NormalCollisions::set_use_area_weighting,
            "If the NormalCollisions should use the convergent formulation.")
        .def_property(
            "collision_set_type", &NormalCollisions::collision_set_type,
            &NormalCollisions::set_collision_set_type,
            R"ipc_Qu8mg5v7(
            The type of collision set to use.

            This can be either:
              - IPC (Implicit Potential Collisions)
              - IMPROVED_MAX_APPROX (Improved Max Approximation)
              - OGC (Offset Geometric Contact)
            )ipc_Qu8mg5v7")
        .def_property(
            "enable_shape_derivatives",
            &NormalCollisions::enable_shape_derivatives,
            &NormalCollisions::set_enable_shape_derivatives,
            "If the NormalCollisions are using the convergent formulation.")
        .def_readwrite("vv_collisions", &NormalCollisions::vv_collisions)
        .def_readwrite("ev_collisions", &NormalCollisions::ev_collisions)
        .def_readwrite("ee_collisions", &NormalCollisions::ee_collisions)
        .def_readwrite("fv_collisions", &NormalCollisions::fv_collisions)
        .def_readwrite("pv_collisions", &NormalCollisions::pv_collisions);
}
