#include <common.hpp>

#include <ipc/broad_phase/brute_force.hpp>

namespace py = pybind11;
using namespace ipc;

void define_brute_force(py::module_& m)
{
    py::class_<BruteForce, BroadPhase>(m, "BruteForce").def(py::init());
}
