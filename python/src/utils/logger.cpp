#include <common.hpp>

#include <ipc/utils/logger.hpp>

namespace py = pybind11;
using namespace ipc;

void define_logger(py::module_& m)
{
    py::enum_<spdlog::level::level_enum>(
        m, "LoggerLevel", "Enumeration of log levels")
        .value("trace", spdlog::level::level_enum::trace, "Trace level")
        .value("debug", spdlog::level::level_enum::debug, "Debug level")
        .value("info", spdlog::level::level_enum::info, "Info level")
        .value("warn", spdlog::level::level_enum::warn, "Warning level")
        .value("error", spdlog::level::level_enum::err, "Error level")
        .value(
            "critical", spdlog::level::level_enum::critical, "Critical level")
        .value("off", spdlog::level::level_enum::off, "Off level")
        .export_values();

    m.def(
        "set_logger_level",
        [](const spdlog::level::level_enum& level) {
            logger().set_level(level);
        },
        "Set log level", py::arg("level"));
}
