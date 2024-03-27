#include <common.hpp>

#include <ipc/ccd/ccd.hpp>

#include <tight_inclusion/ccd.hpp>
#include <tight_inclusion/interval_root_finder.hpp>

namespace py = pybind11;
using namespace ipc;

void define_tight_inclusion_ccd(py::module_& m)
{
    auto m_ti = m.def_submodule(
        "tight_inclusion", "Tight Inclusion CCD method of [Wang et al. 2021].");

    py::enum_<ticcd::CCDRootFindingMethod>(
        m_ti, "CCDRootFindingMethod",
        "Enumeration of implemented root finding methods.")
        .value(
            "DEPTH_FIRST_SEARCH",
            ticcd::CCDRootFindingMethod::DEPTH_FIRST_SEARCH,
            "Depth first search")
        .value(
            "BREADTH_FIRST_SEARCH",
            ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH,
            "Breadth first search")
        .export_values();

    m_ti.def(
        "edge_edge_ccd",
        [](const ticcd::Vector3& ea0_t0, const ticcd::Vector3& ea1_t0,
           const ticcd::Vector3& eb0_t0, const ticcd::Vector3& eb1_t0,
           const ticcd::Vector3& ea0_t1, const ticcd::Vector3& ea1_t1,
           const ticcd::Vector3& eb0_t1, const ticcd::Vector3& eb1_t1,
           const ticcd::Scalar min_distance, const ticcd::Scalar tmax,
           const ticcd::Scalar tolerance, const long max_iterations,
           const ticcd::Array3& filter, bool no_zero_toi,
           const ticcd::CCDRootFindingMethod ccd_method) {
            ticcd::Scalar toi;
            ticcd::Scalar output_tolerance;
            bool r = ticcd::edgeEdgeCCD(
                ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1,
                filter, min_distance, toi, tolerance, tmax, max_iterations,
                output_tolerance, no_zero_toi, ccd_method);
            return std::make_tuple(r, toi, output_tolerance);
        },
        R"ipc_Qu8mg5v7(
        Determine the earliest time of impact between two edges (optionally with a minimum separation).

        Parameters:
            ea0_t0: Starting position of the first vertex of the first edge.
            ea1_t0: Start position of the second vertex of the first edge.
            eb0_t0: Start position of the first vertex of the second edge.
            eb1_t0: Start position of the second vertex of the second edge.
            ea0_t1: End position of the first vertex of the first edge.
            ea1_t1: End position of the second vertex of the first edge.
            eb0_t1: End position of the first vertex of the second edge.
            eb1_t1: End position of the second vertex of the second edge.
            min_distance: Minimum separation distance (default: 0).
            tmax: Upper bound of the time interval [0,tmax] to be checked (0<=tmax<=1).
            tolerance: Solver tolerance (default: 1e-6).
            max_iterations: Maximum number of solver iterations (default: 1e7). If negative, solver will run until convergence.
            filter: Filters calculated using get_numerical_error (default: (-1,-1,-1)). Use (-1,-1,-1) if checking a single query.
            no_zero_toi: Refine further if a zero TOI is produced (assuming not initially in contact).
            ccd_method: Root finding method (default: BREADTH_FIRST_SEARCH).

        Returns:
            Tuple of:
                True if there is a collision, false otherwise,
                the earliest time of collision if collision happens (infinity if no collision occurs), and
                if max_iterations < 0, the solver precision otherwise the input tolerance.
        )ipc_Qu8mg5v7",
        py::arg("ea0_t0"), py::arg("ea1_t0"), py::arg("eb0_t0"),
        py::arg("eb1_t0"), py::arg("ea0_t1"), py::arg("ea1_t1"),
        py::arg("eb0_t1"), py::arg("eb1_t1"), py::arg("min_distance") = 0,
        py::arg("tmax") = 1, py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
        py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
        py::arg("filter") = ticcd::Array3::Constant(-1),
        py::arg("no_zero_toi") = ticcd::DEFAULT_NO_ZERO_TOI,
        py::arg("ccd_method") =
            ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    m_ti.def(
        "point_triangle_ccd",
        [](const ticcd::Vector3& v_t0, const ticcd::Vector3& f0_t0,
           const ticcd::Vector3& f1_t0, const ticcd::Vector3& f2_t0,
           const ticcd::Vector3& v_t1, const ticcd::Vector3& f0_t1,
           const ticcd::Vector3& f1_t1, const ticcd::Vector3& f2_t1,
           const ticcd::Scalar min_distance, const ticcd::Scalar tmax,
           const ticcd::Scalar tolerance, const long max_iterations,
           const ticcd::Array3& filter, bool no_zero_toi,
           const ticcd::CCDRootFindingMethod ccd_method) {
            ticcd::Scalar toi;
            ticcd::Scalar output_tolerance;
            bool r = ticcd::vertexFaceCCD(
                v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1, filter,
                min_distance, toi, tolerance, tmax, max_iterations,
                output_tolerance, no_zero_toi, ccd_method);
            return std::make_tuple(r, toi, output_tolerance);
        },
        R"ipc_Qu8mg5v7(
        Determine the earliest time of impact between a point and triangle (optionally with a minimum separation).

        Parameters:
            v_t0:  Starting position of the vertex.
            f0_t0: Starting position of the first vertex of the face.
            f1_t0: Starting position of the second vertex of the face.
            f2_t0: Starting position of the third vertex of the face.
            v_t1:  Ending position of the vertex.
            f0_t1: Ending position of the first vertex of the face.
            f1_t1: Ending position of the second vertex of the face.
            f2_t1: Ending position of the third vertex of the face.
            min_distance: Minimum separation distance (default: 0).
            tmax: Upper bound of the time interval [0,tmax] to be checked (0<=tmax<=1).
            tolerance: Solver tolerance (default: 1e-6).
            max_iterations: Maximum number of solver iterations (default: 1e7). If negative, solver will run until convergence.
            filter: Filters calculated using get_numerical_error (default: (-1,-1,-1)). Use (-1,-1,-1) if checking a single query.
            no_zero_toi: Refine further if a zero TOI is produced (assuming not initially in contact).
            ccd_method: Root finding method (default: BREADTH_FIRST_SEARCH).

        Returns:
            Tuple of:
                True if there is a collision, false otherwise,
                the earliest time of collision if collision happens (infinity if no collision occurs), and
                if max_iterations < 0, the solver precision otherwise the input tolerance.
        )ipc_Qu8mg5v7",
        py::arg("v_t0"), py::arg("f0_t0"), py::arg("f1_t0"), py::arg("f2_t0"),
        py::arg("v_t1"), py::arg("f0_t1"), py::arg("f1_t1"), py::arg("f2_t1"),
        py::arg("min_distance") = 0, py::arg("tmax") = 1,
        py::arg("tolerance") = DEFAULT_CCD_TOLERANCE,
        py::arg("max_iterations") = DEFAULT_CCD_MAX_ITERATIONS,
        py::arg("filter") = ticcd::Array3::Constant(-1),
        py::arg("no_zero_toi") = ticcd::DEFAULT_NO_ZERO_TOI,
        py::arg("ccd_method") =
            ticcd::CCDRootFindingMethod::BREADTH_FIRST_SEARCH);

    m_ti.def(
        "compute_ccd_filters",
        [](const ticcd::Vector3& min_corner, const ticcd::Vector3& max_corner,
           const bool is_vertex_face, const bool using_minimum_separation) {
            std::vector<ticcd::Vector3> vertices(2);
            vertices[0] = min_corner;
            vertices[1] = max_corner;
            return ticcd::get_numerical_error(
                vertices, is_vertex_face, using_minimum_separation);
        },
        R"ipc_Qu8mg5v7(
        Compute the numerical error filters for the input to the CCD solver.

        Before you run the simulation, you need to conservatively estimate the
        axis-aligned bounding box in which the meshes will be located during the
        whole simulation process.

        Parameters:
            min_corner: Minimum corner of the axis-aligned bounding box of the simulation scene.
            max_corner: Maximum corner of the axis-aligned bounding box of the simulation scene.
            is_vertex_face: True if checking vertex-face collision, false if checking edge-edge collision.
            using_minimum_separation: True if using minimum separation CCD, false otherwise.

        Returns:
            The numerical error filters for the input parameters.
        )ipc_Qu8mg5v7",
        py::arg("min_corner"), py::arg("max_corner"), py::arg("is_vertex_face"),
        py::arg("using_minimum_separation"));
}