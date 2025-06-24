#include <common.hpp>

#include <ipc/ccd/point_static_plane.hpp>

using namespace ipc;

void define_point_static_plane(py::module_& m)
{
    m.def(
        "point_static_plane_ccd",
        [](Eigen::ConstRef<VectorMax3d> p_t0, Eigen::ConstRef<VectorMax3d> p_t1,
           Eigen::ConstRef<VectorMax3d> plane_origin,
           Eigen::ConstRef<VectorMax3d> plane_normal,
           const double conservative_rescaling) {
            double toi;
            bool r = point_static_plane_ccd(
                p_t0, p_t1, plane_origin, plane_normal, toi,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        R"ipc_Qu8mg5v7(
        Compute the time of impact between a point and a static plane in 3D using continuous collision detection.

        Parameters:
            p_t0: The initial position of the point.
            p_t1: The final position of the point.
            plane_origin: The origin of the plane.
            plane_normal: The normal of the plane.
            conservative_rescaling: Conservative rescaling of the time of impact.

        Returns:
            Tuple of:
            True if a collision was detected, false otherwise.
            Output time of impact
        )ipc_Qu8mg5v7",
        "p_t0"_a, "p_t1"_a, "plane_origin"_a, "plane_normal"_a,
        "conservative_rescaling"_a =
            TightInclusionCCD::DEFAULT_CONSERVATIVE_RESCALING);
}
