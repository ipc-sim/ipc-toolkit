#include <common.hpp>

#include <ipc/config.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>
#ifdef IPC_TOOLKIT_WITH_GLTF
#include <ipc/io/write_gltf.hpp>
#include <ipc/io/read_gltf.hpp>
#endif

#include <pybind11/detail/common.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace ipc;
using namespace ipc::rigid;

PYBIND11_MAKE_OPAQUE(std::vector<Pose>)

void define_rigid_bodies(py::module_& m)
{
    // module_local(false): stl_bind types are module-local by default, which
    // would prevent any other extension module embedding the toolkit from
    // returning Poses to Python.
    py::bind_vector<std::vector<Pose>>(m, "Poses", py::module_local(false))
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
                Eigen::ConstRef<VectorMax3d>,
                Eigen::ConstRef<VectorMax3d>>(),
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

    py::class_<RigidBody> rigid_body(m, "RigidBody");

    py::enum_<RigidBody::Type>(
        rigid_body, "Type", "How a body is simulated.")
        .value("STATIC", RigidBody::Type::STATIC, "Does not move.")
        .value(
            "KINEMATIC", RigidBody::Type::KINEMATIC,
            "Moves along a prescribed script/velocity but does not respond "
            "to forces.")
        .value(
            "DYNAMIC", RigidBody::Type::DYNAMIC,
            "Moves and responds to forces.")
        .export_values();

    rigid_body
        .def(
            py::init<
                Eigen::Ref<Eigen::MatrixXd>, Eigen::ConstRef<Eigen::MatrixXi>,
                Eigen::ConstRef<Eigen::MatrixXi>, const double, Pose&,
                const VectorMax6b&>(),
            "vertices"_a, "edges"_a, "faces"_a, "density"_a, "initial_pose"_a,
            "is_dof_fixed"_a = VectorMax6b())
        .def_property_readonly("mass", &RigidBody::mass)
        .def_property_readonly(
            "moment_of_inertia", &RigidBody::moment_of_inertia)
        .def_property_readonly("J", &RigidBody::J)
        .def_property_readonly("R0", &RigidBody::R0)
        .def_property_readonly("external_force", &RigidBody::external_force)
        .def_property(
            "type", &RigidBody::type, &RigidBody::set_type,
            "How this body is simulated (STATIC/KINEMATIC/DYNAMIC).")
        .def_property_readonly(
            "is_dof_fixed", &RigidBody::is_dof_fixed,
            "Per-DOF fixed flags ([position | rotation]; set at "
            "construction).")
        .def(
            "convert_to_static", &RigidBody::convert_to_static,
            "Convert this body to a STATIC body (all DOF fixed). "
            "Kinematic scripting/driving lives on the demo Simulator "
            "(see ipctk.demo.KinematicDriver / Simulator.set_kinematic_driver).");

    py::class_<RigidBodies, std::shared_ptr<RigidBodies>, CollisionMesh>(
        m, "RigidBodies")
        .def(
            py::init(&RigidBodies::build_from_meshes), "rest_positions"_a,
            "edges"_a, "faces"_a, "densities"_a, "initial_poses"_a,
            "convert_planes"_a = false,
            "is_dof_fixed"_a = std::vector<VectorMax6b>())
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

#ifdef IPC_TOOLKIT_WITH_GLTF
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
#endif
}
