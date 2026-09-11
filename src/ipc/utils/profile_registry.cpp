#include "profile_registry.hpp"

#include <ipc/utils/logger.hpp>

#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

namespace ipc {

ProfileRegistry& ProfileRegistry::instance()
{
    static ProfileRegistry registry;
    return registry;
}

void ProfileRegistry::add_time(const std::string& name, double ms)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    auto& s = m_stats[name];
    s.total += ms;
    s.count += 1;
}

void ProfileRegistry::add_value(const std::string& name, double value)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    auto& s = m_stats[name];
    s.total += value;
    s.count += 1;
}

void ProfileRegistry::reset()
{
    std::lock_guard<std::mutex> lock(m_mutex);
    m_stats.clear();
}

void ProfileRegistry::dump_json(const std::string& path) const
{
    // Snapshot the map under the lock, then release before touching disk so
    // the hot path (add_time / add_value) is not blocked on file I/O.
    // std::map for stable alphabetical key order in the JSON output.
    std::map<std::string, Stat> sorted_copy;
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        for (const auto& [name, s] : m_stats) {
            sorted_copy.emplace(name, s);
        }
    }

    std::ostringstream out;
    out << "{\n";
    bool first = true;
    out << std::setprecision(12);
    for (const auto& [name, s] : sorted_copy) {
        if (!first) {
            out << ",\n";
        }
        first = false;
        const double mean =
            s.count > 0 ? s.total / static_cast<double>(s.count) : 0.0;
        out << "  \"" << name << "\": {"
            << "\"total\": " << s.total << ", "
            << "\"count\": " << s.count << ", "
            << "\"mean\": " << mean << "}";
    }
    out << "\n}\n";

    std::ofstream f(path);
    if (!f.is_open()) {
        logger().error(
            "ProfileRegistry::dump_json: failed to open '{}' for writing",
            path);
        return;
    }
    f << out.str();
}

} // namespace ipc
