#include <common.hpp>

#include <ipc/broad_phase/bvh.hpp>

using namespace ipc;

void define_bvh(py::module_& m)
{
    py::class_<BVH, BroadPhase, std::shared_ptr<BVH>>(m, "BVH").def(py::init());
}
