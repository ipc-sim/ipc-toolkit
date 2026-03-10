#include <common.hpp>

#include <ipc/dynamics/rigid/simulator.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#include <ipc/io/write_gltf.hpp>
#include <ipc/io/read_gltf.hpp>

#include <pybind11/detail/common.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace ipc;
using namespace ipc::rigid;

PYBIND11_MAKE_OPAQUE(std::vector<Pose>)

void define_rigid_simulator(py::module_& m)
{
    py::bind_vector<std::vector<Pose>>(m, "Poses")
        .def("__repr__", [](const std::vector<Pose>& poses) {
            std::string repr = "Poses([";
            for (const auto& pose : poses) {
                repr += fmt::format(
                    "Pose(position={}, rotation={}){}", pose.position,
                    pose.rotation, (&pose != &poses.back() ? ", " : ""));
            }
            repr += "])";
            return repr;
        });

    py::class_<Pose>(m, "Pose")
        .def(py::init<>())
        .def(
            py::init<
                Eigen::ConstRef<Eigen::Vector3d>,
                Eigen::ConstRef<Eigen::Vector3d>>(),
            "position"_a, "rotation"_a)
        .def("rotation_matrix", &Pose::rotation_matrix)
        .def("transform_vertices", &Pose::transform_vertices, "vertices"_a)
        .def(py::self * py::self)
        .def_readwrite("position", &Pose::position)
        .def_readwrite("rotation", &Pose::rotation)
        .def("__repr__", [](const Pose& p) {
            return fmt::format(
                "Pose(position={}, rotation={})", p.position, p.rotation);
        });

    py::class_<RigidBody>(m, "RigidBody")
        .def(
            py::init<
                Eigen::Ref<Eigen::MatrixXd>, Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>, const double, Pose&>(),
            "vertices"_a, "edges"_a, "faces"_a, "density"_a, "initial_pose"_a)
        .def_property_readonly("mass", &RigidBody::mass)
        .def_property_readonly(
            "moment_of_inertia", &RigidBody::moment_of_inertia)
        .def_property_readonly("J", &RigidBody::J)
        .def_property_readonly("R0", &RigidBody::R0)
        .def_property_readonly("external_force", &RigidBody::external_force);

    py::class_<RigidBodies, std::shared_ptr<RigidBodies>, CollisionMesh>(
        m, "RigidBodies")
        .def(
            py::init(&RigidBodies::build_from_meshes), "rest_positions"_a,
            "edges"_a, "faces"_a, "densities"_a, "initial_poses"_a,
            "convert_planes"_a = false)
        .def(
            "vertices",
            py::overload_cast<const std::vector<Pose>&>(
                &RigidBodies::vertices, py::const_),
            R"ipc_Qu8mg5v7(
             Compute the vertex positions from the poses of the rigid bodies.

             Parameters:
                 poses: The poses of the rigid bodies.

             Returns:
                 The vertex positions of the rigid bodies (#V × dim).
             )ipc_Qu8mg5v7",
            "poses"_a)
        .def_property_readonly("num_bodies", &RigidBodies::num_bodies)
        .def("__len__", &RigidBodies::num_bodies)
        .def(
            "__getitem__", py::overload_cast<size_t>(&RigidBodies::operator[]),
            "index"_a);

    py::class_<Simulator>(m, "Simulator")
        .def(
            py::init<
                const std::shared_ptr<RigidBodies>, const std::vector<Pose>&,
                const double>(),
            "bodies"_a, "initial_poses"_a, "dt"_a)
        .def("run", &rigid::Simulator::run, "t_end"_a, "callback"_a)
        .def("step", &rigid::Simulator::step)
        .def("reset", &rigid::Simulator::reset)
        .def_property_readonly(
            "pose_history", &rigid::Simulator::pose_history,
            R"ipc_Qu8mg5v7(
             Get the history of poses in the simulation.

             Returns:
                 A list of poses at each time step.
             )ipc_Qu8mg5v7")
        .def_property_readonly("t", &rigid::Simulator::t);

    m.def(
        "write_gltf", &rigid::write_gltf, R"ipc_Qu8mg5v7(
         Write a sequence of rigid body poses to a glTF file.

         Parameters:
             filename: The output glTF filename.
             bodies: The rigid bodies to write.
             poses: A list of poses for each timestep.
             timestep: The time interval between each pose in seconds.
             embed_buffers: Whether to embed the binary buffers in the glTF file.
             write_binary: Whether to write a binary .glb file (true) or a text .gltf file (false).
             prettyPrint: Whether to pretty-print the JSON content.

         Returns:
             True if successful, false otherwise.
         )ipc_Qu8mg5v7",
        "filename"_a, "bodies"_a, "poses"_a, "timestep"_a,
        "embed_buffers"_a = true, "write_binary"_a = true,
        "prettyPrint"_a = true);

    m.def(
        "read_gltf", &rigid::read_gltf, R"ipc_Qu8mg5v7(
         Read a rigid body scene from a glTF file and return a RigidBodies object.

         Parameters:
             filename: The input glTF filename.
             convert_planes: Whether to convert plane primitives in the glTF file to infinite planes in the RigidBodies object. If false, plane primitives will be ignored.

         Returns:
             A pair containing the RigidBodies object and a vector of initial poses for each body.
         )ipc_Qu8mg5v7",
        "filename"_a, "convert_planes"_a = false);
}
