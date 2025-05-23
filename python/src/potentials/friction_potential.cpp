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
            "The smooth friction mollifier parameter :math:`\\epsilon_{v}`.")
        .def_property(
            "material_friction_table", 
            &FrictionPotential::get_material_friction_table, 
            &FrictionPotential::set_material_friction_table,
            "The material friction coefficient lookup table.")
        .def(
            "gradient",
            [](const FrictionPotential& self, 
               const Eigen::MatrixXd& V, 
               const Eigen::VectorXd& u0, 
               const TangentBasis& basis, 
               double dt, double mu, double eps_v,
               double s_mu, double k_mu) {
                return self.gradient(V, u0, basis, dt, mu, eps_v, s_mu, k_mu);
            },
            R"ipc_Qu8mg5v7(
            Compute the gradient of the friction potential with static and kinetic friction.

            Parameters:
                V: Vertices positions
                u0: Displacement
                basis: Tangent basis
                dt: Time step
                mu: Friction coefficient (if s_mu and k_mu not provided)
                eps_v: Smoothing parameter
                s_mu: Static friction coefficient
                k_mu: Kinetic friction coefficient
            )ipc_Qu8mg5v7",
            py::arg("V"), py::arg("u0"), py::arg("basis"), py::arg("dt"),
            py::arg("mu"), py::arg("eps_v"), py::arg("s_mu") = -1.0, py::arg("k_mu") = -1.0);
}
