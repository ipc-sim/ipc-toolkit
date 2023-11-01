#include <common.hpp>

#include <ipc/distance/edge_edge_mollifier.hpp>

namespace py = pybind11;
using namespace ipc;

void define_edge_edge_mollifier(py::module_& m)
{
    m.def(
        "edge_edge_cross_squarednorm", &edge_edge_cross_squarednorm,
        R"ipc_Qu8mg5v7(
        Compute the squared norm of the edge-edge cross product.

        Parameters:
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.

        Returns:
            The squared norm of the edge-edge cross product.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_cross_squarednorm_gradient",
        py::overload_cast<
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&>(
            &edge_edge_cross_squarednorm_gradient),
        R"ipc_Qu8mg5v7(
        Compute the gradient of the squared norm of the edge cross product.

        Parameters:
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.

        Returns:
            The gradient of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_cross_squarednorm_hessian",
        py::overload_cast<
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&>(
            &edge_edge_cross_squarednorm_hessian),
        R"ipc_Qu8mg5v7(
        Compute the hessian of the squared norm of the edge cross product.

        Parameters:
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.

        Returns:
            The hessian of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_mollifier",
        py::overload_cast<const double, const double>(&edge_edge_mollifier),
        R"ipc_Qu8mg5v7(
        Mollifier function for edge-edge distance.

        Parameters:
            x: Squared norm of the edge-edge cross product.
            eps_x: Mollifier activation threshold.

        Returns:
            The mollifier coefficient to premultiply the edge-edge distance.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier_gradient",
        py::overload_cast<const double, const double>(
            &edge_edge_mollifier_gradient),
        R"ipc_Qu8mg5v7(
        The gradient of the mollifier function for edge-edge distance.

        Parameters:
            x: Squared norm of the edge-edge cross product.
            eps_x: Mollifier activation threshold.

        Returns:
            The gradient of the mollifier function for edge-edge distance wrt x.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier_derivative_wrt_eps_x",
        &edge_edge_mollifier_derivative_wrt_eps_x,
        R"ipc_Qu8mg5v7(
        The derivative of the mollifier function for edge-edge distance wrt eps_x.

        Parameters:
            x: Squared norm of the edge-edge cross product.
            eps_x: Mollifier activation threshold.

        Returns:
            The derivative of the mollifier function for edge-edge distance wrt eps_x.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier_hessian",
        py::overload_cast<const double, const double>(
            &edge_edge_mollifier_hessian),
        R"ipc_Qu8mg5v7(
        The hessian of the mollifier function for edge-edge distance.

        Parameters:
            x: Squared norm of the edge-edge cross product.
            eps_x: Mollifier activation threshold.

        Returns:
            The hessian of the mollifier function for edge-edge distance wrt x.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier_gradient_derivative_wrt_eps_x",
        &edge_edge_mollifier_gradient_derivative_wrt_eps_x,
        R"ipc_Qu8mg5v7(
        The derivative of the gradient of the mollifier function for edge-edge distance wrt eps_x.

        Parameters:
            x: Squared norm of the edge-edge cross product.
            eps_x: Mollifier activation threshold.

        Returns:
            The derivative of the gradient of the mollifier function for edge-edge distance wrt eps_x.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier",
        py::overload_cast<
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&, const double>(
            &edge_edge_mollifier),
        R"ipc_Qu8mg5v7(
        Compute a mollifier for the edge-edge distance.

        This helps smooth the non-smoothness at close to parallel edges.

        Parameters:
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.
            eps_x: Mollifier activation threshold.

        Returns:
            The mollifier coefficient to premultiply the edge-edge distance.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier_gradient",
        py::overload_cast<
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&, const double>(
            &edge_edge_mollifier_gradient),
        R"ipc_Qu8mg5v7(
        Compute the gradient of the mollifier for the edge-edge distance.

        Parameters:
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.
            eps_x: Mollifier activation threshold.

        Returns:
            The gradient of the mollifier.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier_hessian",
        py::overload_cast<
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&,
            const Eigen::Ref<const Eigen::Vector3d>&, const double>(
            &edge_edge_mollifier_hessian),
        R"ipc_Qu8mg5v7(
        Compute the hessian of the mollifier for the edge-edge distance.

        Parameters:
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.
            eps_x: Mollifier activation threshold.

        Returns:
            The hessian of the mollifier.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier_gradient_wrt_x",
        &edge_edge_mollifier_gradient_wrt_x,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the mollifier for the edge-edge distance wrt rest positions.

        Parameters:
            ea0_rest: The rest position of the first vertex of the first edge.
            ea1_rest: The rest position of the second vertex of the first edge.
            eb0_rest: The rest position of the first vertex of the second edge.
            eb1_rest: The rest position of the second vertex of the second edge.
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.

        Returns:
            The derivative of the mollifier wrt rest positions.
        )ipc_Qu8mg5v7",
        py::arg("ea0_rest"), py::arg("ea1_rest"), py::arg("eb0_rest"),
        py::arg("eb1_rest"), py::arg("ea0"), py::arg("ea1"), py::arg("eb0"),
        py::arg("eb1"));

    m.def(
        "edge_edge_mollifier_gradient_jacobian_wrt_x",
        &edge_edge_mollifier_gradient_jacobian_wrt_x,
        R"ipc_Qu8mg5v7(
        Compute the jacobian of the edge-edge distance mollifier's gradient wrt rest positions.

        Note:
            This is not the hessian of the mollifier wrt rest positions, but the jacobian wrt rest positions of the mollifier's gradient wrt positions.

        Parameters:
            ea0_rest: The rest position of the first vertex of the first edge.
            ea1_rest: The rest position of the second vertex of the first edge.
            eb0_rest: The rest position of the first vertex of the second edge.
            eb1_rest: The rest position of the second vertex of the second edge.
            ea0: The first vertex of the first edge.
            ea1: The second vertex of the first edge.
            eb0: The first vertex of the second edge.
            eb1: The second vertex of the second edge.

        Returns:
            The jacobian of the mollifier's gradient wrt rest positions.
        )ipc_Qu8mg5v7",
        py::arg("ea0_rest"), py::arg("ea1_rest"), py::arg("eb0_rest"),
        py::arg("eb1_rest"), py::arg("ea0"), py::arg("ea1"), py::arg("eb0"),
        py::arg("eb1"));

    m.def(
        "edge_edge_mollifier_threshold", &edge_edge_mollifier_threshold,
        R"ipc_Qu8mg5v7(
        Compute the threshold of the mollifier edge-edge distance.

        This values is computed based on the edges at rest length.

        Parameters:
            ea0_rest: The rest position of the first vertex of the first edge.
            ea1_rest: The rest position of the second vertex of the first edge.
            eb0_rest: The rest position of the first vertex of the second edge.
            eb1_rest: The rest position of the second vertex of the second edge.

        Returns:
            Threshold for edge-edge mollification.
        )ipc_Qu8mg5v7",
        py::arg("ea0_rest"), py::arg("ea1_rest"), py::arg("eb0_rest"),
        py::arg("eb1_rest"));

    m.def(
        "edge_edge_mollifier_threshold_gradient",
        &edge_edge_mollifier_threshold_gradient,
        R"ipc_Qu8mg5v7(
        Compute the gradient of the threshold of the mollifier edge-edge distance.

        This values is computed based on the edges at rest length.

        Parameters:
            ea0_rest: The rest position of the first vertex of the first edge.
            ea1_rest: The rest position of the second vertex of the first edge.
            eb0_rest: The rest position of the first vertex of the second edge.
            eb1_rest: The rest position of the second vertex of the second edge.

        Returns:
            Gradient of the threshold for edge-edge mollification.
        )ipc_Qu8mg5v7",
        py::arg("ea0_rest"), py::arg("ea1_rest"), py::arg("eb0_rest"),
        py::arg("eb1_rest"));
}
