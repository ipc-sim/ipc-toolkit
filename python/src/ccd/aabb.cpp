#include <common.hpp>

#include <ipc/ccd/aabb.hpp>

using namespace ipc;

void define_ccd_aabb(py::module_& m)
{
    m.def(
        "point_edge_aabb_cd", &point_edge_aabb_cd, "p"_a, "e0"_a, "e1"_a,
        "dist"_a);

    m.def(
        "edge_edge_aabb_cd", &edge_edge_aabb_cd, "ea0"_a, "ea1"_a, "eb0"_a,
        "eb1"_a, "dist"_a);

    m.def(
        "point_triangle_aabb_cd", &point_triangle_aabb_cd, "p"_a, "t0"_a,
        "t1"_a, "t2"_a, "dist"_a);

    m.def(
        "edge_triangle_aabb_cd", &edge_triangle_aabb_cd, "e0"_a, "e1"_a, "t0"_a,
        "t1"_a, "t2"_a, "dist"_a);

    m.def(
        "point_edge_aabb_ccd", &point_edge_aabb_ccd, "p_t0"_a, "e0_t0"_a,
        "e1_t0"_a, "p_t1"_a, "e0_t1"_a, "e1_t1"_a, "dist"_a);

    m.def(
        "edge_edge_aabb_ccd", &edge_edge_aabb_ccd, "ea0_t0"_a, "ea1_t0"_a,
        "eb0_t0"_a, "eb1_t0"_a, "ea0_t1"_a, "ea1_t1"_a, "eb0_t1"_a, "eb1_t1"_a,
        "dist"_a);

    m.def(
        "point_triangle_aabb_ccd", &point_triangle_aabb_ccd, "p_t0"_a,
        "t0_t0"_a, "t1_t0"_a, "t2_t0"_a, "p_t1"_a, "t0_t1"_a, "t1_t1"_a,
        "t2_t1"_a, "dist"_a);
}
