#include <common.hpp>

#include <ipc/broad_phase/bvh.hpp>

namespace py = pybind11;
using namespace ipc;

void define_bvh(py::module_& m)
{
    py::class_<BVH, BroadPhase>(m, "BVH").def(py::init());
}
