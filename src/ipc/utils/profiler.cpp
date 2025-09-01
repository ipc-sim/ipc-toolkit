#include "profiler.hpp"

#ifdef IPC_TOOLKIT_WITH_PROFILER

namespace ipc {

Profiler& profiler()
{
    static Profiler instance;
    return instance;
}

} // namespace ipc

#endif