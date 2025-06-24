#include <common.hpp>

#include <ipc/ccd/inexact_point_edge.hpp>

using namespace ipc;

void define_inexact_point_edge(py::module_& m)
{
    m.def(
        "inexact_point_edge_ccd_2D",
        [](Eigen::ConstRef<Eigen::Vector2d> p_t0,
           Eigen::ConstRef<Eigen::Vector2d> e0_t0,
           Eigen::ConstRef<Eigen::Vector2d> e1_t0,
           Eigen::ConstRef<Eigen::Vector2d> p_t1,
           Eigen::ConstRef<Eigen::Vector2d> e0_t1,
           Eigen::ConstRef<Eigen::Vector2d> e1_t1,
           const double conservative_rescaling) {
            double toi;
            bool r = inexact_point_edge_ccd_2D(
                p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Inexact continuous collision detection between a point and an edge in 2D.

        Parameters:
            p_t0: Initial position of the point
            e0_t0: Initial position of the first endpoint of the edge
            e1_t0: Initial position of the second endpoint of the edge
            p_t1: Final position of the point
            e0_t1: Final position of the first endpoint of the edge
            e1_t1: Final position of the second endpoint of the edge
            conservative_rescaling: Conservative rescaling of the time of impact

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            Output time of impact
        )ipc_Qu8mg5v7",
        "p_t0"_a, "e0_t0"_a, "e1_t0"_a, "p_t1"_a, "e0_t1"_a, "e1_t1"_a,
        "conservative_rescaling"_a);
}
