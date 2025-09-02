#include "profiler.hpp"

#ifdef IPC_TOOLKIT_WITH_PROFILER

#include <fstream>

namespace ipc {

Profiler& profiler()
{
    static Profiler instance;
    return instance;
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
    os << "Parent,Name,Time (ms),Count\n";

    if (m_data.empty()) {
        os << std::flush;
        return;
    }

    // Print the profiler data in CSV format using a breadth-first traversal
    const nlohmann::json::json_pointer root;
    std::queue<nlohmann::json::json_pointer> queue;
    queue.push(root);

    while (!queue.empty()) {
        const auto ptr = queue.front();
        const auto parent = ptr.parent_pointer();
        queue.pop();

        assert(m_data.contains(ptr));
        const auto& data = ptr == root ? m_data : m_data.at(ptr);
        if (ptr != root) {
            os << fmt::format(
                "{},{},{:.6g},{:d}\n", parent == root ? "" : parent.back(),
                ptr.back(), data.at("time_ms").get<double>(),
                data.at("count").get<size_t>());
        }

        // Traverse child scopes
        for (const auto& [key, val] : data.items()) {
            if (val.is_object()) {
                queue.push(ptr / key);
            }
        }
    }

    os << std::flush;
}

} // namespace ipc

#endif