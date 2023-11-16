#include <common.hpp>

#include <ipc/ccd/point_static_plane.hpp>

namespace py = pybind11;
using namespace ipc;

void define_point_static_plane(py::module_& m)
{
    m.def(
        "point_static_plane_ccd",
        [](const VectorMax3d& p_t0, const VectorMax3d& p_t1,
           const VectorMax3d& plane_origin, const VectorMax3d& plane_normal,
           const double conservative_rescaling =
               DEFAULT_CCD_CONSERVATIVE_RESCALING) {
            double toi;
            bool r = point_static_plane_ccd(
                p_t0, p_t1, plane_origin, plane_normal, toi,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        py::arg("p_t0"), py::arg("p_t1"), py::arg("plane_origin"),
        py::arg("plane_normal"),
        py::arg("conservative_rescaling") = DEFAULT_CCD_CONSERVATIVE_RESCALING);
}
