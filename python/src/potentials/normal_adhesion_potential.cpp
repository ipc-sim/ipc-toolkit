#include <common.hpp>

#include <ipc/potentials/normal_adhesion_potential.hpp>

namespace py = pybind11;
using namespace ipc;

void define_normal_adhesion_potential(py::module_& m)
{
    py::class_<NormalAdhesionPotential, NormalPotential>(
        m, "NormalAdhesionPotential")
        .def(
            py::init<const double, const double, const double, const double>(),
            py::arg("dhat_p"), py::arg("dhat_a"), py::arg("Y"),
            py::arg("eps_c"))
        .def_readwrite(
            "dhat_p", &NormalAdhesionPotential::dhat_p,
            "The distance of largest adhesion force (:math:`\\hat{d}_{p}`) (:math:`0 < \\hat{d}_{p} < \\hat{d}_{a}`).")
        .def_readwrite(
            "dhat_a", &NormalAdhesionPotential::dhat_a,
            "The adhesion activation distance (:math:`\\hat{d}_{a}`.")
        .def_readwrite(
            "Y", &NormalAdhesionPotential::Y,
            "The Young's modulus (:math:`Y`)).")
        .def_readwrite(
            "eps_c", &NormalAdhesionPotential::eps_c,
            "The critical strain (:math:`\\varepsilon_{c}`)).");
}
