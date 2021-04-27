// clang-format off
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>
// clang-format on

#include <tbb/global_control.h>
#include <tbb/task_scheduler_init.h>

#include <ipc/barrier/barrier.hpp>
#include <ipc/collision_constraint.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/distance_type.hpp>
#include <ipc/distance/edge_edge_mollifier.hpp>
#include <ipc/distance/edge_edge.hpp>
#include <ipc/distance/line_line.hpp>
#include <ipc/distance/point_line.hpp>
#include <ipc/distance/point_plane.hpp>
#include <ipc/distance/point_point.hpp>
#include <ipc/distance/point_triangle.hpp>
#include <ipc/utils/logger.hpp>

#include "utils.hpp"

namespace py = pybind11;
using namespace ipc;

// static tbb::global_control thread_limiter(
//     tbb::global_control::max_allowed_parallelism,
//     tbb::task_scheduler_init::default_num_threads());
//
// void set_num_threads(int nthreads)
// {
//     if (nthreads <= 0) {
//         nthreads = tbb::task_scheduler_init::default_num_threads();
//     } else if (nthreads > tbb::task_scheduler_init::default_num_threads()) {
//         logger().warn(
//             "Attempting to use more threads than available ({:d} > "
//             "{:d})!",
//             nthreads, tbb::task_scheduler_init::default_num_threads());
//         nthreads = tbb::task_scheduler_init::default_num_threads();
//     }
//     thread_limiter = tbb::global_control(
//         tbb::global_control::max_allowed_parallelism, nthreads);
// }

PYBIND11_MODULE(ipctk, m)
{
    m.doc() = "IPC Toolkit";

    ///////////////////////////////////////////////////////////////////////////
    // barrier/barrier
    m.def(
        "barrier", &barrier<double>,
        R"ipc_Qu8mg5v7(
        Function that grows to infinity as d approaches 0 from the right.

        Parameters
        ----------
        d : The distance
        dhat : Activation distance of the barrier

        Returns
        -------
        The value of the barrier function at d.

        Notes
        -----
        .. math:: b(d) = -(d-\hat{d})^2\ln(\frac{d}{\hat{d}})
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));

    m.def(
        "barrier_gradient", &barrier_gradient,
        R"ipc_Qu8mg5v7(
        Derivative of the barrier function.

        Parameters
        ----------
        d : The distance
        dhat : Activation distance of the barrier

        Returns
        -------
        The derivative of the barrier wrt d.

        Notes
        -----
        .. math:: b(d) = (\hat{d}-d)(2\ln(\frac{d}{\hat{d}}) - \frac{\hat{d}}{d}) + 1)
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));

    m.def(
        "barrier_hessian", &barrier_hessian,
        R"ipc_Qu8mg5v7(
        Second derivative of the barrier function.

        Parameters
        ----------
        d : The distance
        dhat : Activation distance of the barrier

        Returns
        -------
        The second derivative of the barrier wrt d.

        Notes
        -----
        .. math:: b(d) = (\frac{\hat{d}}{d} + 2)\frac{\hat{d}}{d} - 2\ln(\frac{d}{\hat{d}}) - 3
        )ipc_Qu8mg5v7",
        py::arg("d"), py::arg("dhat"));

    ///////////////////////////////////////////////////////////////////////////
    // collision_constraint
    py::class_<VertexVertexCandidate>(m, "VertexVertexCandidate")
        .def(py::init<long, long>())
        .def_readwrite("vertex0_index", &VertexVertexConstraint::vertex0_index)
        .def_readwrite("vertex1_index", &VertexVertexConstraint::vertex1_index);
    py::class_<EdgeVertexCandidate>(m, "EdgeVertexCandidate")
        .def(py::init<long, long>())
        .def(
            "__str__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format(
                    "[{:d}, {:d}]", ev.edge_index, ev.vertex_index);
            })
        .def(
            "__repr__",
            [](const EdgeVertexCandidate& ev) {
                return fmt::format(
                    "EdgeVertexCandidate({:d}, {:d})", ev.edge_index,
                    ev.vertex_index);
            })
        .def_readwrite("edge_index", &EdgeVertexCandidate::edge_index)
        .def_readwrite("vertex_index", &EdgeVertexCandidate::vertex_index);
    py::class_<EdgeEdgeCandidate>(m, "EdgeEdgeCandidate")
        .def(py::init<long, long>())
        .def_readwrite("edge0_index", &EdgeEdgeCandidate::edge0_index)
        .def_readwrite("edge1_index", &EdgeEdgeCandidate::edge1_index);
    py::class_<EdgeFaceCandidate>(m, "EdgeFaceCandidate")
        .def(py::init<long, long>())
        .def_readwrite("edge_index", &EdgeFaceCandidate::edge_index)
        .def_readwrite("face_index", &EdgeFaceCandidate::face_index);
    py::class_<FaceVertexCandidate>(m, "FaceVertexCandidate")
        .def(py::init<long, long>())
        .def_readwrite("face_index", &FaceVertexCandidate::face_index)
        .def_readwrite("vertex_index", &FaceVertexCandidate::vertex_index);
    py::class_<Candidates>(m, "Candidates")
        .def(py::init())
        .def("size", &Candidates::size)
        .def("clear", &Candidates::clear)
        .def_readwrite("ev_candidates", &Candidates::ev_candidates)
        .def_readwrite("ee_candidates", &Candidates::ee_candidates)
        .def_readwrite("fv_candidates", &Candidates::fv_candidates);

    py::class_<VertexVertexConstraint, VertexVertexCandidate>(
        m, "VertexVertexConstraint")
        .def(py::init<long, long>())
        // .def(py::init<const VertexVertexCandidate&>())
        .def(
            "vertex_indices", &VertexVertexConstraint::vertex_indices,
            R"ipc_Qu8mg5v7(
            Get the indices of the vertices

            Parameters
            ----------
            E : Edge matrix of mesh
            F : Face matrix of mesh

            Returns
            -------
            List of vertex indices
            )ipc_Qu8mg5v7",
            py::arg("E"), py::arg("F"))
        .def(
            "compute_potential", &VertexVertexConstraint::compute_potential,
            py::arg("V"), py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_gradient",
            &VertexVertexConstraint::compute_potential_gradient, py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"))
        .def(
            "compute_potential_hessian",
            &VertexVertexConstraint::compute_potential_hessian, py::arg("V"),
            py::arg("E"), py::arg("F"), py::arg("dhat"),
            py::arg("project_to_psd"))
        .def_readwrite("multiplicity", &VertexVertexConstraint::multiplicity);

    ///////////////////////////////////////////////////////////////////////////
    // distance/distance_type

    py::enum_<PointEdgeDistanceType>(m, "PointEdgeDistanceType")
        .value("P_E0", PointEdgeDistanceType::P_E0)
        .value("P_E1", PointEdgeDistanceType::P_E1)
        .value("P_E", PointEdgeDistanceType::P_E)
        .export_values();

    py::enum_<PointTriangleDistanceType>(m, "PointTriangleDistanceType")
        .value("P_T0", PointTriangleDistanceType::P_T0)
        .value("P_T1", PointTriangleDistanceType::P_T1)
        .value("P_T2", PointTriangleDistanceType::P_T2)
        .value("P_E0", PointTriangleDistanceType::P_E0)
        .value("P_E1", PointTriangleDistanceType::P_E1)
        .value("P_E2", PointTriangleDistanceType::P_E2)
        .value("P_T", PointTriangleDistanceType::P_T)
        .export_values();

    py::enum_<EdgeEdgeDistanceType>(m, "EdgeEdgeDistanceType")
        .value("EA0_EB0", EdgeEdgeDistanceType::EA0_EB0)
        .value("EA0_EB1", EdgeEdgeDistanceType::EA0_EB1)
        .value("EA1_EB0", EdgeEdgeDistanceType::EA1_EB0)
        .value("EA1_EB1", EdgeEdgeDistanceType::EA1_EB1)
        .value("EA_EB0", EdgeEdgeDistanceType::EA_EB0)
        .value("EA_EB1", EdgeEdgeDistanceType::EA_EB1)
        .value("EA0_EB", EdgeEdgeDistanceType::EA0_EB)
        .value("EA1_EB", EdgeEdgeDistanceType::EA1_EB)
        .value("EA_EB", EdgeEdgeDistanceType::EA_EB)
        .export_values();

    m.def(
        "point_edge_distance_type",
        [](const Eigen::VectorX3d& p, const Eigen::VectorX3d& e0,
           const Eigen::VectorX3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_edge_distance_type(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Determine the closest pair between a point and edge.

        Parameters
        ----------
        p  : The point
        e0 : The first vertex of the edge
        e1 : The second vertex of the edge

        Returns
        -------
        The distance type of the point-edge pair.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_triangle_distance_type",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            return point_triangle_distance_type(p, t0, t1, t2);
        },
        R"ipc_Qu8mg5v7(
        Determine the closest pair between a point and triangle.

        Parameters
        ----------
        p  : The point
        t0 : The first vertex of the triangle
        t1 : The second vertex of the triangle
        t2 : The third vertex of the triangle

        Returns
        -------
        The distance type of the point-triangle pair.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"));

    m.def(
        "edge_edge_distance_type",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return edge_edge_distance_type(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Determine the closest pair between two edges.

        Parameters
        ----------
        ea0 : The first vertex of the first edge
        ea1 : The second vertex of the first edge
        eb0 : The first vertex of the second edge
        eb1 : The second vertex of the second edge

        Returns
        -------
        The distance type of the edge-edge pair.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    ///////////////////////////////////////////////////////////////////////////
    // distance/edge_edge_mollifier

    m.def(
        "edge_edge_cross_squarednorm",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return edge_edge_cross_squarednorm(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Compute the squared norm of the edge-edge cross product.

        Parameters
        ----------
        ea0 : The first vertex of the first edge
        ea1 : The second vertex of the first edge
        eb0 : The first vertex of the second edge
        eb1 : The second vertex of the second edge

        Returns
        -------
        The squared norm of the edge-edge cross product.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_cross_squarednorm_gradient",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            Eigen::Vector<double, 12> grad;
            edge_edge_cross_squarednorm_gradient(ea0, ea1, eb0, eb1, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the squared norm of the edge cross product.

        Parameters
        ----------
        ea0 : The first vertex of the first edge
        ea1 : The second vertex of the first edge
        eb0 : The first vertex of the second edge
        eb1 : The second vertex of the second edge

        Returns
        -------
        The gradient of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_cross_squarednorm_hessian",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            Eigen::Matrix<double, 12, 12> hess;
            edge_edge_cross_squarednorm_hessian(ea0, ea1, eb0, eb1, hess);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the squared norm of the edge cross product.

        Parameters
        ----------
        ea0 : The first vertex of the first edge
        ea1 : The second vertex of the first edge
        eb0 : The first vertex of the second edge
        eb1 : The second vertex of the second edge

        Returns
        -------
        The hessian of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "edge_edge_mollifier", &edge_edge_mollifier<double>,
        R"ipc_Qu8mg5v7(
        Mollifier function for edge-edge distance.

        Parameters
        ----------
        x : Squared norm of the edge-edge cross product
        eps_x : Mollifier activation threshold.

        Returns
        -------
        The mollifier coefficient to premultiply the edge-edge distance.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier", &edge_edge_mollifier_gradient<double>,
        R"ipc_Qu8mg5v7(
        The gradient of the mollifier function for edge-edge distance.

        Parameters
        ----------
        x : Squared norm of the edge-edge cross product
        eps_x : Mollifier activation threshold.

        Returns
        -------
        The gradient of the mollifier function for edge-edge distance wrt x.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier", &edge_edge_mollifier_hessian<double>,
        R"ipc_Qu8mg5v7(
        The hessian of the mollifier function for edge-edge distance.

        Parameters
        ----------
        x : Squared norm of the edge-edge cross product
        eps_x : Mollifier activation threshold.

        Returns
        -------
        The hessian of the mollifier function for edge-edge distance wrt x.
        )ipc_Qu8mg5v7",
        py::arg("x"), py::arg("eps_x"));

    m.def(
        "edge_edge_mollifier_threshold",
        [](const Eigen::Vector3d& ea0_rest, const Eigen::Vector3d& ea1_rest,
           const Eigen::Vector3d& eb0_rest, const Eigen::Vector3d& eb1_rest) {
            return edge_edge_mollifier_threshold(
                ea0_rest, ea1_rest, eb0_rest, eb1_rest);
        },
        R"ipc_Qu8mg5v7(
        Compute the threshold of the mollifier edge-edge distance.

        This values is computed based on the edges at rest length.

        Parameters
        ----------
        x : Squared norm of the edge-edge cross product
        eps_x : Mollifier activation threshold.

        Returns
        -------
        Threshold for edge-edge mollification.
        )ipc_Qu8mg5v7",
        py::arg("ea0_rest"), py::arg("ea1_rest"), py::arg("eb0_rest"),
        py::arg("eb1_rest"));

    ///////////////////////////////////////////////////////////////////////////
    // distance/edge_edge

    m.def(
        "edge_edge_distance",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1,
           const EdgeEdgeDistanceType* dtype) {
            if (dtype == nullptr) {
                return edge_edge_distance(ea0, ea1, eb0, eb1);
            } else {
                return edge_edge_distance(ea0, ea1, eb0, eb1, *dtype);
            }
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a two lines segments in 3D.

        Parameters
        ----------
        ea0 : The first vertex of the first edge
        ea1 : The second vertex of the first edge
        eb0 : The first vertex of the second edge
        eb1 : The second vertex of the second edge
        dtype : (Optional) The edge-edge distance type to compute

        Returns
        -------
        The distance between the two edges.

        See also
        --------
        edge_edge_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("dtype") = py::none());

    m.def(
        "edge_edge_distance_gradient",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1,
           const EdgeEdgeDistanceType* dtype) {
            Eigen::Vector<double, 12> grad;
            if (dtype == nullptr) {
                edge_edge_distance_gradient(ea0, ea1, eb0, eb1, grad);
            } else {
                edge_edge_distance_gradient(ea0, ea1, eb0, eb1, *dtype, grad);
            }
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a two lines segments.

        Parameters
        ----------
        ea0 : The first vertex of the first edge
        ea1 : The second vertex of the first edge
        eb0 : The first vertex of the second edge
        eb1 : The second vertex of the second edge
        dtype : (Optional) The point edge distance type to compute

        Returns
        -------
        The gradient of the distance wrt ea0, ea1, eb0, and eb1.

        See also
        --------
        edge_edge_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("dtype") = py::none());

    m.def(
        "edge_edge_distance_hessian",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1,
           const EdgeEdgeDistanceType* dtype) {
            Eigen::Matrix<double, 12, 12> hess;
            if (dtype == nullptr) {
                edge_edge_distance_hessian(ea0, ea1, eb0, eb1, hess);
            } else {
                edge_edge_distance_hessian(ea0, ea1, eb0, eb1, *dtype, hess);
            }
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a two lines segments.

        Parameters
        ----------
        ea0 : The first vertex of the first edge
        ea1 : The second vertex of the first edge
        eb0 : The first vertex of the second edge
        eb1 : The second vertex of the second edge
        dtype : (Optional) The point edge distance type to compute

        Returns
        -------
        The hessian of the distance wrt ea0, ea1, eb0, and eb1.

        See also
        --------
        edge_edge_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"),
        py::arg("dtype") = py::none());

    ///////////////////////////////////////////////////////////////////////////
    // distance/line_line

    m.def(
        "line_line_distance",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            return line_line_distance(ea0, ea1, eb0, eb1);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a two infinite lines in 3D.

        Parameters
        ----------
        ea0 : The first vertex of the edge defining the first line
        ea1 : The second vertex of the edge defining the first line
        eb0 : The first vertex of the edge defining the second line
        eb1 : The second vertex of the edge defining the second line

        Returns
        -------
        The distance between the two lines.

        Notes
        -----
        The distance is actually squared distance.

        Warning
        -------
        If the lines are parallel this function returns a distance of zero.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "line_line_distance_gradient",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            Eigen::Vector<double, 12> grad;
            line_line_distance_gradient(ea0, ea1, eb0, eb1, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a two lines in 3D.

        Parameters
        ----------
        ea0 : The first vertex of the edge defining the first line
        ea1 : The second vertex of the edge defining the first line
        eb0 : The first vertex of the edge defining the second line
        eb1 : The second vertex of the edge defining the second line

        Returns
        -------
        The gradient of the distance wrt ea0, ea1, eb0, and eb1.

        Notes
        -----
        The distance is actually squared distance.

        Warning
        -------
        If the lines are parallel this function returns a distance of zero.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    m.def(
        "line_line_distance_hessian",
        [](const Eigen::Vector3d& ea0, const Eigen::Vector3d& ea1,
           const Eigen::Vector3d& eb0, const Eigen::Vector3d& eb1) {
            Eigen::Matrix<double, 12, 12> hess;
            line_line_distance_hessian(ea0, ea1, eb0, eb1, hess);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a two lines in 3D.

        Parameters
        ----------
        ea0 : The first vertex of the edge defining the first line
        ea1 : The second vertex of the edge defining the first line
        eb0 : The first vertex of the edge defining the second line
        eb1 : The second vertex of the edge defining the second line

        Returns
        -------
        The hessian of the distance wrt ea0, ea1, eb0, and eb1.

        Notes
        -----
        The distance is actually squared distance.

        Warning
        -------
        If the lines are parallel this function returns a distance of zero.
        )ipc_Qu8mg5v7",
        py::arg("ea0"), py::arg("ea1"), py::arg("eb0"), py::arg("eb1"));

    ///////////////////////////////////////////////////////////////////////////
    // distance/point_edge

    m.def(
        "point_edge_distance",
        [](const Eigen::VectorX3d& p, const Eigen::VectorX3d& e0,
           const Eigen::VectorX3d& e1, const PointEdgeDistanceType* dtype) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            if (dtype == nullptr) {
                return point_edge_distance(p, e0, e1);
            } else {
                return point_edge_distance(p, e0, e1, *dtype);
            }
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a point and edge in 2D or 3D.

        Parameters
        ----------
        p : The point
        e0 : The first vertex of the edge
        e1 : The second vertex of the edge
        dtype : (Optional) The point edge distance type to compute

        Returns
        -------
        The distance between the point and edge

        See also
        --------
        point_edge_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = py::none());

    m.def(
        "point_edge_distance_gradient",
        [](const Eigen::VectorX3d& p, const Eigen::VectorX3d& e0,
           const Eigen::VectorX3d& e1, const PointEdgeDistanceType* dtype) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            Eigen::VectorX9d grad;
            if (dtype == nullptr) {
                point_edge_distance_gradient(p, e0, e1, grad);
            } else {
                point_edge_distance_gradient(p, e0, e1, *dtype, grad);
            }
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a point and edge.

        Parameters
        ----------
        p : The point
        e0 : The first vertex of the edge
        e1 : The second vertex of the edge
        dtype : (Optional) The point edge distance type to compute

        Returns
        -------
        The gradient of the distance wrt p, e0, and e1.

        See also
        --------
        point_edge_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = py::none());

    m.def(
        "point_edge_distance_hessian",
        [](const Eigen::VectorX3d& p, const Eigen::VectorX3d& e0,
           const Eigen::VectorX3d& e1, const PointEdgeDistanceType* dtype) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            Eigen::MatrixXX9d hess;
            if (dtype == nullptr) {
                point_edge_distance_hessian(p, e0, e1, hess);
            } else {
                point_edge_distance_hessian(p, e0, e1, *dtype, hess);
            }
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and edge.

        Parameters
        ----------
        p : The point
        e0 : The first vertex of the edge
        e1 : The second vertex of the edge
        dtype : (Optional) The point edge distance type to compute

        Returns
        -------
        The hessian of the distance wrt p, e0, and e1.

        See also
        --------
        point_edge_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"),
        py::arg("dtype") = py::none());

    ///////////////////////////////////////////////////////////////////////////
    // distance/point_line

    m.def(
        "point_line_distance",
        [](const Eigen::VectorX3d& p, const Eigen::VectorX3d& e0,
           const Eigen::VectorX3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            return point_line_distance(p, e0, e1);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a point and line in 2D or 3D.

        Parameters
        ----------
        p : The point
        e0 : The first vertex of the edge defining the line
        e1 : The second vertex of the edge defining the line

        Returns
        -------
        The distance between the point and line.

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_line_distance_gradient",
        [](const Eigen::VectorX3d& p, const Eigen::VectorX3d& e0,
           const Eigen::VectorX3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            Eigen::VectorX9d grad;
            point_line_distance_gradient(p, e0, e1, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a point and line.

        Parameters
        ----------
        p : The point
        e0 : The first vertex of the edge defining the line.
        e1 : The second vertex of the edge defining the line.

        Returns
        -------
        The gradient of the distance wrt p, e0, and e1.

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    m.def(
        "point_line_distance_hessian",
        [](const Eigen::VectorX3d& p, const Eigen::VectorX3d& e0,
           const Eigen::VectorX3d& e1) {
            assert_2D_or_3D_vector(p, "p");
            assert_2D_or_3D_vector(e0, "e0");
            assert_2D_or_3D_vector(e1, "e1");
            Eigen::MatrixXX9d hess;
            point_line_distance_hessian(p, e0, e1, hess);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and line.

        Parameters
        ----------
        p : The point
        e0 : The first vertex of the edge defining the line
        e1 : The second vertex of the edge defining the line

        Returns
        -------
        The hessian of the distance wrt p, e0, and e1.

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("e0"), py::arg("e1"));

    ///////////////////////////////////////////////////////////////////////////
    // distance/point_plane

    m.def(
        "point_plane_distance",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            return point_plane_distance(p, t0, t1, t2);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a point and a plane.

        Parameters
        ----------
        p  : The point.
        t0 : The first vertex of the triangle.
        t1 : The second vertex of the triangle.
        t2 : The third vertex of the triangle.

        Returns
        -------
        The distance between the point and plane.

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t1"));

    m.def(
        "point_plane_distance_gradient",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            Eigen::Vector<double, 12> grad;
            point_plane_distance_gradient(p, t0, t1, t2, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a point and a plane.

        Parameters
        ----------
        p  : The point
        t0 : The first vertex of the triangle
        t1 : The second vertex of the triangle
        t2 : The third vertex of the triangle

        Returns
        -------
        The gradient of the distance wrt p, t0, t1, and t2.

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t1"));

    m.def(
        "point_plane_distance_hessian",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
            Eigen::Matrix<double, 12, 12> hess;
            point_plane_distance_hessian(p, t0, t1, t2, hess);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and a plane.

        Parameters
        ----------
        p  : The point
        t0 : The first vertex of the triangle
        t1 : The second vertex of the triangle
        t2 : The third vertex of the triangle

        Returns
        -------
        The hessian of the distance wrt p, t0, t1, and t2.

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t1"));

    ///////////////////////////////////////////////////////////////////////////
    // distance/point_point

    m.def(
        "point_point_distance",
        [](const Eigen::VectorX3d& p0, const Eigen::VectorX3d& p1) {
            assert_2D_or_3D_vector(p0, "p0");
            assert_2D_or_3D_vector(p1, "p1");
            return point_point_distance(p0, p1);
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between two points.

        Parameters
        ----------
        p0 : The first point
        p1 : The second point

        Returns
        -------
        The distance between p0 and p1

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p0"), py::arg("p1"));

    m.def(
        "point_point_distance_gradient",
        [](const Eigen::VectorX3d& p0, const Eigen::VectorX3d& p1) {
            assert_2D_or_3D_vector(p0, "p0");
            assert_2D_or_3D_vector(p1, "p1");
            Eigen::VectorX6<double> grad;
            point_point_distance_gradient(p0, p1, grad);
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between two points.

        Parameters
        ----------
        p0 : The first point
        p1 : The second point

        Returns
        -------
        The gradient of the distance wrt p0 and p1.

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p0"), py::arg("p1"));

    m.def(
        "point_point_distance_hessian",
        [](const Eigen::VectorX3d& p0, const Eigen::VectorX3d& p1) {
            assert_2D_or_3D_vector(p0, "p0");
            assert_2D_or_3D_vector(p1, "p1");
            Eigen::MatrixXX6<double> hess;
            point_point_distance_hessian(p0, p1, hess);
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a point and point.

        Parameters
        ----------
        p0 : The first point
        p1 : The second point

        Returns
        -------
        The hessian of the distance wrt p0 and p1.

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p0"), py::arg("p1"));

    ///////////////////////////////////////////////////////////////////////////
    // distance/point_triangle

    m.def(
        "point_triangle_distance",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2,
           const PointTriangleDistanceType* dtype) {
            if (dtype == nullptr) {
                return point_triangle_distance(p, t0, t1, t2);
            } else {
                return point_triangle_distance(p, t0, t1, t2, *dtype);
            }
        },
        R"ipc_Qu8mg5v7(
        Compute the distance between a two lines segments in 3D.

        Parameters
        ----------
        p  : The point.
        t0 : The first vertex of the triangle.
        t1 : The second vertex of the triangle.
        t2 : The third vertex of the triangle.
        dtype : (Optional) The point-triangle distance type to compute

        Returns
        -------
        The distance between the point and triangle.

        See also
        --------
        point_triangle_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"),
        py::arg("dtype") = py::none());

    m.def(
        "point_triangle_distance_gradient",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2,
           const PointTriangleDistanceType* dtype) {
            Eigen::Vector<double, 12> grad;
            if (dtype == nullptr) {
                point_triangle_distance_gradient(p, t0, t1, t2, grad);
            } else {
                point_triangle_distance_gradient(p, t0, t1, t2, *dtype, grad);
            }
            return grad;
        },
        R"ipc_Qu8mg5v7(
        Compute the gradient of the distance between a two lines segments.

        Parameters
        ----------
        p  : The point.
        t0 : The first vertex of the triangle.
        t1 : The second vertex of the triangle.
        t2 : The third vertex of the triangle.
        dtype : (Optional) The point-triangle distance type to compute

        Returns
        -------
        The gradient of the distance wrt ea0, ea1, eb0, and eb1.

        See also
        --------
        point_triangle_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"),
        py::arg("dtype") = py::none());

    m.def(
        "point_triangle_distance_hessian",
        [](const Eigen::Vector3d& p, const Eigen::Vector3d& t0,
           const Eigen::Vector3d& t1, const Eigen::Vector3d& t2,
           const PointTriangleDistanceType* dtype) {
            Eigen::Matrix<double, 12, 12> hess;
            if (dtype == nullptr) {
                point_triangle_distance_hessian(p, t0, t1, t2, hess);
            } else {
                point_triangle_distance_hessian(p, t0, t1, t2, *dtype, hess);
            }
            return hess;
        },
        R"ipc_Qu8mg5v7(
        Compute the hessian of the distance between a two lines segments.

        Parameters
        ----------
        p  : The point.
        t0 : The first vertex of the triangle.
        t1 : The second vertex of the triangle.
        t2 : The third vertex of the triangle.
        dtype : (Optional) The point-triangle distance type to compute

        Returns
        -------
        The hessian of the distance wrt ea0, ea1, eb0, and eb1.

        See also
        --------
        point_triangle_distance_type

        Notes
        -----
        The distance is actually squared distance.
        )ipc_Qu8mg5v7",
        py::arg("p"), py::arg("t0"), py::arg("t1"), py::arg("t2"),
        py::arg("dtype") = py::none());

    ///////////////////////////////////////////////////////////////////////////
    // utils/logger

    py::enum_<spdlog::level::level_enum>(m, "LoggerLevel")
        .value("trace", spdlog::level::level_enum::trace)
        .value("debug", spdlog::level::level_enum::debug)
        .value("info", spdlog::level::level_enum::info)
        .value("warn", spdlog::level::level_enum::warn)
        .value("error", spdlog::level::level_enum::err)
        .value("critical", spdlog::level::level_enum::critical)
        .value("off", spdlog::level::level_enum::off)
        .export_values();

    m.def(
        "set_logger_level",
        [](const spdlog::level::level_enum& level) {
            logger().set_level(level);
        },
        "Set log level", py::arg("level"));

    // m.def(
    //     "set_num_threads", &set_num_threads, "maximum number of threads to
    //     use", py::arg("nthreads"));
}
