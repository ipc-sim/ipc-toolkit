#include <common.hpp>

#include <ipc/implicits/plane.hpp>

namespace py = pybind11;
using namespace ipc;

void define_plane_implicit(py::module_& m)
{
    m.def(
        "construct_point_plane_collisions",
        [](const Eigen::MatrixXd& points, const Eigen::MatrixXd& plane_origins,
           const Eigen::MatrixXd& plane_normals, const double dhat,
           const double dmin = 0) {
            std::vector<PlaneVertexCollision> pv_collisions;
            construct_point_plane_collisions(
                points, plane_origins, plane_normals, dhat, pv_collisions,
                dmin);
            return pv_collisions;
        },
        R"ipc_Qu8mg5v7(
        Construct a set of point-plane distance collisions used to compute

        Note:
            The given pv_collisions will be cleared.

        the barrier potential.

        Parameters:
            points: Points as rows of a matrix.
            plane_origins: Plane origins as rows of a matrix.
            plane_normals: Plane normals as rows of a matrix.
            dhat: The activation distance of the barrier.
            dmin: Minimum distance.

        Returns:
            The constructed set of collisions.
        )ipc_Qu8mg5v7",
        py::arg("points"), py::arg("plane_origins"), py::arg("plane_normals"),
        py::arg("dhat"), py::arg("dmin") = 0);

    m.def(
        "construct_point_plane_collisions",
        [](const Eigen::MatrixXd& points, const Eigen::MatrixXd& plane_origins,
           const Eigen::MatrixXd& plane_normals, const double dhat,
           const double dmin,
           const std::function<bool(size_t, size_t)>& can_collide) {
            std::vector<PlaneVertexCollision> pv_collisions;
            construct_point_plane_collisions(
                points, plane_origins, plane_normals, dhat, pv_collisions, dmin,
                can_collide);
            return pv_collisions;
        },
        R"ipc_Qu8mg5v7(
        Construct a set of point-plane distance collisions used to compute

        Note:
            The given pv_collisions will be cleared.

        the barrier potential.

        Parameters:
            points: Points as rows of a matrix.
            plane_origins: Plane origins as rows of a matrix.
            plane_normals: Plane normals as rows of a matrix.
            dhat: The activation distance of the barrier.
            dmin: Minimum distance.
            can_collide: A function that takes a vertex ID (row numbers in points) and a plane ID (row number in plane_origins) then returns true if the vertex can collide with the plane. By default all points can collide with all planes.

        Returns:
            The constructed set of collisions.
        )ipc_Qu8mg5v7",
        py::arg("points"), py::arg("plane_origins"), py::arg("plane_normals"),
        py::arg("dhat"), py::arg("dmin"), py::arg("can_collide"));

    m.def(
        "is_step_point_plane_collision_free",
        [](const Eigen::MatrixXd& points_t0, const Eigen::MatrixXd& points_t1,
           const Eigen::MatrixXd& plane_origins,
           const Eigen::MatrixXd& plane_normals) {
            return is_step_point_plane_collision_free(
                points_t0, points_t1, plane_origins, plane_normals);
        },
        R"ipc_Qu8mg5v7(
        Determine if the step is collision free.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            points_t0: Points at start as rows of a matrix.
            points_t1: Points at end as rows of a matrix.
            plane_origins: Plane origins as rows of a matrix.
            plane_normals: Plane normals as rows of a matrix.

        Returns:
            True if <b>any</b> collisions occur.
        )ipc_Qu8mg5v7",
        py::arg("points_t0"), py::arg("points_t1"), py::arg("plane_origins"),
        py::arg("plane_normals"));

    m.def(
        "is_step_point_plane_collision_free",
        [](const Eigen::MatrixXd& points_t0, const Eigen::MatrixXd& points_t1,
           const Eigen::MatrixXd& plane_origins,
           const Eigen::MatrixXd& plane_normals,
           const std::function<bool(size_t, size_t)>& can_collide) {
            return is_step_point_plane_collision_free(
                points_t0, points_t1, plane_origins, plane_normals,
                can_collide);
        },
        R"ipc_Qu8mg5v7(
        Determine if the step is collision free.

        Note:
            Assumes the trajectory is linear.

        Parameters:
            points_t0: Points at start as rows of a matrix.
            points_t1: Points at end as rows of a matrix.
            plane_origins: Plane origins as rows of a matrix.
            plane_normals: Plane normals as rows of a matrix.
            can_collide: A function that takes a vertex ID (row numbers in points) and a plane ID (row number in plane_origins) then returns true if the vertex can collide with the plane. By default all points can collide with all planes.

        Returns:
            True if <b>any</b> collisions occur.
        )ipc_Qu8mg5v7",
        py::arg("points_t0"), py::arg("points_t1"), py::arg("plane_origins"),
        py::arg("plane_normals"), py::arg("can_collide"));

    m.def(
        "compute_point_plane_collision_free_stepsize",
        [](const Eigen::MatrixXd& points_t0, const Eigen::MatrixXd& points_t1,
           const Eigen::MatrixXd& plane_origins,
           const Eigen::MatrixXd& plane_normals) {
            return compute_point_plane_collision_free_stepsize(
                points_t0, points_t1, plane_origins, plane_normals);
        },
        R"ipc_Qu8mg5v7(
        Computes a maximal step size that is collision free.

        Notes:
            Assumes points_t0 is intersection free.
            Assumes the trajectory is linear.
            A value of 1.0 if a full step and 0.0 is no step.

        Parameters:
            points_t0: Points at start as rows of a matrix.
            points_t1: Points at end as rows of a matrix.
            plane_origins: Plane origins as rows of a matrix.
            plane_normals: Plane normals as rows of a matrix.
            can_collide: A function that takes a vertex ID (row numbers in points) and a plane ID (row number in plane_origins) then returns true if the vertex can collide with the plane. By default all points can collide with all planes.

        Returns:
            A step-size $\in [0, 1]$ that is collision free.
        )ipc_Qu8mg5v7",
        py::arg("points_t0"), py::arg("points_t1"), py::arg("plane_origins"),
        py::arg("plane_normals"));

    m.def(
        "compute_point_plane_collision_free_stepsize",
        [](const Eigen::MatrixXd& points_t0, const Eigen::MatrixXd& points_t1,
           const Eigen::MatrixXd& plane_origins,
           const Eigen::MatrixXd& plane_normals,
           const std::function<bool(size_t, size_t)>& can_collide) {
            return compute_point_plane_collision_free_stepsize(
                points_t0, points_t1, plane_origins, plane_normals,
                can_collide);
        },
        R"ipc_Qu8mg5v7(
        Computes a maximal step size that is collision free.

        Notes:
            Assumes points_t0 is intersection free.
            Assumes the trajectory is linear.
            A value of 1.0 if a full step and 0.0 is no step.

        Parameters:
            points_t0: Points at start as rows of a matrix.
            points_t1: Points at end as rows of a matrix.
            plane_origins: Plane origins as rows of a matrix.
            plane_normals: Plane normals as rows of a matrix.
            can_collide: A function that takes a vertex ID (row numbers in points) and a plane ID (row number in plane_origins) then returns true if the vertex can collide with the plane. By default all points can collide with all planes.

        Returns:
            A step-size $\in [0, 1]$ that is collision free.
        )ipc_Qu8mg5v7",
        py::arg("points_t0"), py::arg("points_t1"), py::arg("plane_origins"),
        py::arg("plane_normals"), py::arg("can_collide"));
}
