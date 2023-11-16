#include <common.hpp>

#include <ipc/ccd/inexact_point_edge.hpp>

namespace py = pybind11;
using namespace ipc;

void define_inexact_point_edge(py::module_& m)
{
    m.def(
        "inexact_point_edge_ccd_2D",
        [](const Eigen::Vector2d& p_t0, const Eigen::Vector2d& e0_t0,
           const Eigen::Vector2d& e1_t0, const Eigen::Vector2d& p_t1,
           const Eigen::Vector2d& e0_t1, const Eigen::Vector2d& e1_t1,
           const double conservative_rescaling) {
            double toi;
            bool r = inexact_point_edge_ccd_2D(
                p_t0, e0_t0, e1_t0, p_t1, e0_t1, e1_t1, toi,
                conservative_rescaling);
            return std::make_tuple(r, toi);
        },
        py::arg("p_t0"), py::arg("e0_t0"), py::arg("e1_t0"), py::arg("p_t1"),
        py::arg("e0_t1"), py::arg("e1_t1"), py::arg("conservative_rescaling"));
}
