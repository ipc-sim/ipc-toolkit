#include <common.hpp>

#include <ipc/dynamics/rigid/simulator.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

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
            py::init<const Eigen::Vector3d&, const Eigen::Vector3d&>(),
            py::arg("position"), py::arg("rotation"))
        .def("rotation_matrix", &Pose::rotation_matrix)
        .def(
            "transform_vertices", &Pose::transform_vertices,
            py::arg("vertices"))
        .def(py::self * py::self)
        .def_readwrite("position", &Pose::position)
        .def_readwrite("rotation", &Pose::rotation)
        .def("__repr__", [](const Pose& p) {
            return fmt::format(
                "Pose(position={}, rotation={})", p.position, p.rotation);
        });

    py::class_<RigidBodies, std::shared_ptr<RigidBodies>, CollisionMesh>(
        m, "RigidBodies")
        .def(
            py::init(&RigidBodies::build_from_meshes),
            py::arg("rest_positions"), py::arg("edges"), py::arg("faces"),
            py::arg("densities"), py::arg("initial_poses"))
        .def(
            "vertices", &RigidBodies::vertices,
            R"ipc_Qu8mg5v7(
             Compute the vertex positions from the poses of the rigid bodies.

             Parameters:
                 poses: The poses of the rigid bodies.

             Returns:
                 The vertex positions of the rigid bodies (#V × dim).
             )ipc_Qu8mg5v7",
            py::arg("poses"))
        .def_property_readonly("num_bodies", &RigidBodies::num_bodies)
        .def("__getitem__", &RigidBodies::operator[], py::arg("index"));

    py::class_<Simulator>(m, "Simulator")
        .def(
            py::init<
                const std::shared_ptr<RigidBodies>, const std::vector<Pose>&,
                const double>(),
            py::arg("bodies"), py::arg("initial_poses"), py::arg("dt"))
        .def(
            "run", &rigid::Simulator::run, py::arg("t_end"),
            py::arg("callback"))
        .def("step", &rigid::Simulator::step)
        .def("reset", &rigid::Simulator::reset)
        .def_property_readonly(
            "poses_history", &rigid::Simulator::poses_history,
            R"ipc_Qu8mg5v7(
             Get the history of poses in the simulation.

             Returns:
                 A list of poses at each time step.
             )ipc_Qu8mg5v7")
        .def_property_readonly("t", &rigid::Simulator::t);
}
