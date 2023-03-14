#include <common.hpp>

#include <ipc/utils/logger.hpp>

namespace py = pybind11;
using namespace ipc;

void define_logger(py::module_& m)
{
    py::enum_<spdlog::level::level_enum>(m, "LoggerLevel")
        .value("trace", spdlog::level::level_enum::trace)
        .value("debug", spdlog::level::level_enum::debug)
        .value("info", spdlog::level::level_enum::info)
        .value("warn", spdlog::level::level_enum::warn)
        .value("error", spdlog::level::level_enum::err)
        .value("critical", spdlog::level::level_enum::critical)
        .value("off", spdlog::level::level_enum::off)
        .export_values();

    m.def(
        "set_logger_level",
        [](const spdlog::level::level_enum& level) {
            logger().set_level(level);
        },
        "Set log level", py::arg("level"));
}
