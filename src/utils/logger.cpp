#ifdef IPC_TOOLKIT_WITH_LOGGER
#include <iostream>
#include <memory>
#include <mutex>

#include <ipc/utils/logger.hpp>

#include <spdlog/details/registry.h>
#include <spdlog/details/thread_pool.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

namespace ipc {

std::shared_ptr<spdlog::async_logger> Logger::_logger;

// See https://github.com/gabime/spdlog#asynchronous-logger-with-multi-sinks
void Logger::init(std::vector<spdlog::sink_ptr>& sinks)
{
    auto l = spdlog::get("ipc_toolkit");
    bool had_ipc_toolkit = l != nullptr;
    if (had_ipc_toolkit) {
        spdlog::drop("ipc_toolkit");
    }

    if (spdlog::thread_pool() == nullptr) {
        spdlog::init_thread_pool(8192, 1);
    }
    Logger::_logger = std::make_shared<spdlog::async_logger>(
        "ipc_toolkit", sinks.begin(), sinks.end(), spdlog::thread_pool(),
        spdlog::async_overflow_policy::block);
    spdlog::register_logger(_logger);

    if (had_ipc_toolkit) {
        logger().warn("Removed another ipc toolkit logger");
    }
}

void Logger::init(bool use_cout, const std::string& filename, bool truncate)
{
    std::vector<spdlog::sink_ptr> sinks;
    if (use_cout) {
        sinks.emplace_back(
            std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
    }
    if (!filename.empty()) {
        sinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(
            filename, truncate));
    }

    init(sinks);
}

void Logger::init(std::ostream& os)
{
    std::vector<spdlog::sink_ptr> sinks;
    sinks.emplace_back(
        std::make_shared<spdlog::sinks::ostream_sink_mt>(os, false));

    init(sinks);
}

} // namespace ipc
#endif
