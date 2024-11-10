#include <common.hpp>

#include <ipc/ccd/check_initial_distance.hpp>

namespace py = pybind11;
using namespace ipc;

void define_check_initial_distance(py::module_& m)
{
    m.def(
        "check_initial_distance",
        [](const double initial_distance, const double min_distance) {
            double toi;
            bool r =
                check_initial_distance(initial_distance, min_distance, toi);
            return std::make_tuple(r, toi);
        },
        py::arg("initial_distance"), py::arg("min_distance"));
}
