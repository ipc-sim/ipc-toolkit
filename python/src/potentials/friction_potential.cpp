#include <common.hpp>

#include <ipc/potentials/friction_potential.hpp>

namespace py = pybind11;
using namespace ipc;

void define_friction_potential(py::module_& m)
{
    py::class_<FrictionPotential, TangentialPotential>(m, "FrictionPotential")
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a friction potential.

            Parameters:
                eps_v: The smooth friction mollifier parameter :math:`\\epsilon_{v}`.
            )ipc_Qu8mg5v7",
            py::arg("eps_v"))
        .def_property(
            "eps_v", &FrictionPotential::eps_v, &FrictionPotential::set_eps_v,
            "The smooth friction mollifier parameter :math:`\\epsilon_{v}`.");
}
