#include <common.hpp>

#include <ipc/ccd/additive_ccd.hpp>

namespace py = pybind11;
using namespace ipc;

void define_additive_ccd(py::module_& m)
{
    py::class_<AdditiveCCD, NarrowPhaseCCD>(m, "AdditiveCCD")
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a new AdditiveCCD object.

            Parameters:
                conservative_rescaling: The conservative rescaling of the time of impact.
            )ipc_Qu8mg5v7",
            py::arg("conservative_rescaling") =
                AdditiveCCD::DEFAULT_CONSERVATIVE_RESCALING)
        .def_static(
            "additive_ccd",
            [](VectorMax12d x, const VectorMax12d& dx,
               const std::function<double(const VectorMax12d&)>&
                   distance_squared,
               const double max_disp_mag, const double min_distance,
               const double tmax, const double conservative_rescaling) {
                double toi;
                bool r = AdditiveCCD::additive_ccd(
                    x, dx, distance_squared, max_disp_mag, toi, min_distance,
                    tmax, conservative_rescaling);
                return std::make_tuple(r, toi);
            },
            R"ipc_Qu8mg5v7(
            Computes the time of impact between two objects using additive continuous collision detection.

            Parameters:
                distance_squared: A function that computes the squared distance between the two objects at a given time.
                min_distance: The minimum distance between the objects.
                tmax: The maximum time to check for collisions.
                conservative_rescaling: The amount to rescale the objects by to ensure conservative advancement.

            Returns:
                Tuple of:
                True if a collision was detected, false otherwise.
                The time of impact between the two objects.
            )ipc_Qu8mg5v7",
            py::arg("x"), py::arg("dx"), py::arg("distance_squared"),
            py::arg("max_disp_mag"), py::arg("min_distance") = 0.0,
            py::arg("tmax") = 1.0,
            py::arg("conservative_rescaling") =
                AdditiveCCD::DEFAULT_CONSERVATIVE_RESCALING)
        .def_readonly_static(
            "DEFAULT_CONSERVATIVE_RESCALING",
            &AdditiveCCD::DEFAULT_CONSERVATIVE_RESCALING,
            "The default conservative rescaling value used to avoid taking steps exactly to impact. Value choosen to based on [Li et al. 2021].")
        .def_readwrite(
            "conservative_rescaling", &AdditiveCCD::conservative_rescaling,
            "The conservative rescaling value used to avoid taking steps exactly to impact.");
}
