#include "profiler.hpp"

#ifdef IPC_TOOLKIT_WITH_PROFILER

#include <fstream>
#include <queue>

namespace ipc {

Profiler::Profiler() { }

Profiler& profiler()
{
    static Profiler instance;
    return instance;
}

void Profiler::clear() { m_data.clear(); }

void Profiler::start(const std::string& name)
{
    current_scope.push_back(name);
    if (!m_data.contains(current_scope)) {
        m_data[current_scope] = {
            { "time_ms", 0 },
            { "count", 0 },
        };
    }
}

void Profiler::stop(const double time_ms)
{
    const static std::string log_fmt_text = fmt::format(
        "[{}] {{}} {{:.6f}} ms",
        fmt::format(fmt::fg(fmt::terminal_color::magenta), "timing"));

    logger().trace(
        fmt::runtime(log_fmt_text), current_scope.to_string(), time_ms);

    assert(m_data.contains(current_scope));
    assert(m_data.at(current_scope).contains("time_ms"));
    assert(m_data.at(current_scope).contains("count"));
    m_data[current_scope]["time_ms"] =
        m_data[current_scope]["time_ms"].get<double>() + time_ms;
    m_data[current_scope]["count"] =
        m_data[current_scope]["count"].get<size_t>() + 1;
    current_scope.pop_back();
}

void Profiler::reset()
{
    m_data.clear();
    // reset the calling thread's scope
    current_scope = nlohmann::json::json_pointer(); // root
}

void Profiler::print() const
{
    logger().info(
        "[{}] profiler: {}",
        fmt::format(fmt::fg(fmt::terminal_color::magenta), "timing"),
        m_data.dump(2));
}

void Profiler::write_csv(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        logger().error("Failed to open file: {}", filename);
        return;
    }
    write_csv(file);
    file.close();
}

void Profiler::write_csv(std::ostream& os) const
{
    os << "Id,Parent,Name,Time (ms),Count\n";

    if (m_data.empty()) {
        os << std::flush;
        return;
    }

    // Print the profiler data in CSV format using a breadth-first traversal
    const nlohmann::json::json_pointer root;
    // parent id, pointer
    std::queue<std::pair<int, nlohmann::json::json_pointer>> queue;
    queue.push(std::make_pair(-1, root));
    int id = -1;

    while (!queue.empty()) {
        const auto [parent_id, ptr] = queue.front();
        queue.pop();

        assert(m_data.contains(ptr));
        const auto& data = ptr == root ? m_data : m_data.at(ptr);
        if (ptr != root) {
            os << fmt::format(
                "{:d},{},{},{:.6g},{:d}\n", id,
                parent_id == -1 ? "" : std::to_string(parent_id), ptr.back(),
                data.at("time_ms").get<double>(),
                data.at("count").get<size_t>());
        }

        // Traverse child scopes
        for (const auto& [key, val] : data.items()) {
            if (val.is_object()) {
                queue.push(std::make_pair(id, ptr / key));
            }
        }

        ++id;
    }

    os << std::flush;
}

} // namespace ipc

#endif