#pragma once

#include <ipc/broad_phase/hash_grid.hpp>

namespace ipc {

inline std::shared_ptr<BroadPhase> make_default_broad_phase()
{
    return std::make_shared<HashGrid>();
}

} // namespace ipc