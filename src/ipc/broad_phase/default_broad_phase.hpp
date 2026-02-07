#pragma once

#include <ipc/broad_phase/lbvh.hpp>

#include <memory>

namespace ipc {

inline std::unique_ptr<BroadPhase> make_default_broad_phase()
{
    return std::make_unique<LBVH>();
}

} // namespace ipc