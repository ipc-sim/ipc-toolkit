#include <common.hpp>

#include <ipc/potentials/tangential_adhesion_potential.hpp>

using namespace ipc;

void define_tangential_adhesion_potential(py::module_& m)
{
    py::class_<TangentialAdhesionPotential, TangentialPotential>(
        m, "TangentialAdhesionPotential")
        .def(
            py::init<const double>(),
            R"ipc_Qu8mg5v7(
            Construct a tangential adhesion potential.

            Parameters:
                eps_a: The tangential adhesion mollifier parameter :math:`\epsilon_a`.
            )ipc_Qu8mg5v7",
            "eps_a"_a)
        .def_property(
            "eps_a", &TangentialAdhesionPotential::eps_a,
            &TangentialAdhesionPotential::set_eps_a,
            "Get the tangential adhesion mollifier parameter :math:`\epsilon_a`.");
}
