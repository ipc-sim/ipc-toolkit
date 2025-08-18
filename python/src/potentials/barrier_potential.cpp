#include <common.hpp>

#include <ipc/potentials/barrier_potential.hpp>

using namespace ipc;

void define_barrier_potential(py::module_& m)
{
    py::class_<BarrierPotential, NormalPotential>(m, "BarrierPotential")
        .def(
            py::init<const double, const bool>(),
            R"ipc_Qu8mg5v7(
            Construct a barrier potential.

            Parameters:
                dhat: The activation distance of the barrier.
            )ipc_Qu8mg5v7",
            "dhat"_a, "use_physical_barrier"_a = false)
        .def(
            py::init<
                const std::shared_ptr<Barrier>, const double, const bool>(),
            R"ipc_Qu8mg5v7(
            Construct a barrier potential.

            Parameters:
                barrier: The barrier function.
                dhat: The activation distance of the barrier.
            )ipc_Qu8mg5v7",
            "barrier"_a, "dhat"_a, "use_physical_barrier"_a = false)
        .def_property(
            "dhat", &BarrierPotential::dhat, &BarrierPotential::set_dhat,
            "Barrier activation distance.")
        .def_property(
            "barrier",
            py::cpp_function(
                &BarrierPotential::barrier, py::return_value_policy::reference),
            &BarrierPotential::set_barrier,
            "Barrier function used to compute the potential.");
}
