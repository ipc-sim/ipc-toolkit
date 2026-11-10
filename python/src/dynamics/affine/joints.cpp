#include <common.hpp>

#include <ipc/dynamics/affine/joints.hpp>
#include <ipc/dynamics/affine/pose.hpp>
#include <ipc/dynamics/rigid/rigid_bodies.hpp>

#include <pybind11/stl_bind.h>

namespace py = pybind11;
using namespace ipc;
using namespace ipc::affine;

// Match the opaque declaration in dynamics/rigid/bodies.cpp so the bound
// "Poses" class is used consistently across translation units.
PYBIND11_MAKE_OPAQUE(std::vector<rigid::Pose>)

void define_affine_joints(py::module_& m)
{
    py::class_<affine::Pose>(m, "AffinePose")
        .def(py::init<>())
        .def_readwrite("position", &affine::Pose::position)
        .def_readwrite("rotation", &affine::Pose::rotation)
        .def("rotation_vector", &affine::Pose::rotation_vector)
        .def(
            "transform_vertices", &affine::Pose::transform_vertices,
            "vertices"_a)
        .def("__repr__", [](const affine::Pose& p) {
            return fmt::format(
                "AffinePose(position={}, rotation={})", p.position,
                p.rotation.reshaped().transpose());
        });

    py::class_<JointConstraints, std::shared_ptr<JointConstraints>>(
        m, "JointConstraints", R"ipc_Qu8mg5v7(
         Linear equality joint constraints for affine bodies [Chen et al. 2022].

         Anchor points are given in world space at the initial configuration
         and converted to material coordinates using the initial poses.
         )ipc_Qu8mg5v7")
        .def(
            py::init<
                const std::shared_ptr<const rigid::RigidBodies>&,
                const std::vector<rigid::Pose>&>(),
            "bodies"_a, "initial_poses"_a)
        .def(
            "add_point_connection", &JointConstraints::add_point_connection,
            "Glue a point of body_a to a point of body_b.", "body_a"_a,
            "body_b"_a, "world_anchor"_a)
        .def(
            "add_fixed_point", &JointConstraints::add_fixed_point,
            "Fix a material point of a body at its initial world position.",
            "body"_a, "world_anchor"_a)
        .def(
            "add_hinge", &JointConstraints::add_hinge,
            "Hinge between two bodies about the axis through two points.",
            "body_a"_a, "body_b"_a, "world_axis_p0"_a, "world_axis_p1"_a)
        .def(
            "add_sliding_plane", &JointConstraints::add_sliding_plane,
            "Constrain a material point of a body to a plane.", "body"_a,
            "world_anchor"_a, "normal"_a)
        .def(
            "add_fixed_body", &JointConstraints::add_fixed_body,
            "Fix all 12 DOFs of a body at their initial values.", "body"_a)
        .def_property_readonly(
            "num_constraints", &JointConstraints::num_constraints)
        .def(
            "residual", &JointConstraints::residual,
            "Evaluate the constraint residual Cx - s.", "x"_a);
}
