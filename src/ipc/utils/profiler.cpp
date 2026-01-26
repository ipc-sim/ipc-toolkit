#include "profiler.hpp"

#ifdef IPC_TOOLKIT_WITH_PROFILER

#include <algorithm>
#include <fstream>
#include <queue>

namespace ipc {

Profiler::Profiler() { m_main_thread_id = std::this_thread::get_id(); }

Profiler& profiler()
{
    static Profiler instance;
    return instance;
}

void Profiler::clear() { m_thread_data.clear(); }

void Profiler::start(const std::string& name)
{
    auto& [m_data, current_scope, scope_path, sync_version] =
        m_thread_data.local();

    if (std::this_thread::get_id() == m_main_thread_id) {
        // --- Main Thread Logic ---
        if (!m_main_thread_data) {
            m_main_thread_data = &m_thread_data.local();
        }

        // SeqLock: Increment to Odd (Write in progress)
        size_t v = m_scope_version.load(std::memory_order_relaxed);
        m_scope_version.store(v + 1, std::memory_order_release);

        scope_path.push_back(name);
        current_scope.push_back(name);

        // SeqLock: Increment to Even (Write finished)
        m_scope_version.store(v + 2, std::memory_order_release);
    } else {
        // --- Worker Thread Logic ---

        size_t global_v = m_scope_version.load(std::memory_order_acquire);

        // If versions mismatch, we must re-sync our scope with the main thread
        if (sync_version != global_v) {

            std::vector<std::string> snapshot_path;
            size_t v_snapshot = 0;

            // SeqLock Read Loop
            while (true) {
                size_t v1 = m_scope_version.load(std::memory_order_acquire);
                if (v1 % 2 != 0) { // Writing?
                    std::this_thread::yield();
                    continue;
                }

                if (m_main_thread_data) {
                    snapshot_path = m_main_thread_data->scope_path;
                } else {
                    snapshot_path.clear();
                }

                std::atomic_thread_fence(std::memory_order_acquire);
                size_t v2 = m_scope_version.load(std::memory_order_acquire);

                if (v1 == v2) {
                    v_snapshot = v1;
                    break;
                }
            }

            // --- Deduplication / Pruning Strategy ---
            // If the main thread is already deep inside the parallel loop, its
            // stack might look like ["Block 3", "Block 4", "Block 5"]. If we
            // are starting "Block 4", we want to inherit ["Block 3"]. Strategy:
            // Search for 'name' in the snapshot from back to front. If found,
            // prune everything after (and including) it.

            auto it =
                std::find(snapshot_path.rbegin(), snapshot_path.rend(), name);
            if (it != snapshot_path.rend()) {
                // Found 'name' in the main stack.
                // The 'base' scope is everything before this occurrence.
                // rbegin() is the last element. rend() is before first.
                // distance from rbegin is how many elements to pop + 1.
                // Or simply: resize to the index of the element.

                // Convert reverse iterator to index
                // forward iterator: it.base() returns iterator to element AFTER
                // the one found so we want to keep everything up to (it.base()
                // - 1) - 1 ?? Simpler: The index of the element found is:
                size_t found_idx =
                    std::distance(begin(snapshot_path), it.base()) - 1;

                // We want to keep [0 ... found_idx - 1]
                snapshot_path.resize(found_idx);
            }

            // Apply snapshot to local state
            scope_path = snapshot_path;
            current_scope = nlohmann::json::json_pointer();
            for (const auto& s : scope_path) {
                current_scope.push_back(s);
            }

            sync_version = v_snapshot;
        }

        // Finally push the new block for this thread
        scope_path.push_back(name);
        current_scope.push_back(name);
    }

    // Initialize data container if needed
    if (!m_data.contains(current_scope)) {
        m_data[current_scope] = {
            { "time_ms", 0.0 },
            { "count", 0 },
        };
    }
}

void Profiler::stop(const double time_ms)
{
    const static std::string log_fmt_text = fmt::format(
        "[{}] {{}} {{:.6f}} ms",
        fmt::format(fmt::fg(fmt::terminal_color::magenta), "timing"));

    auto& [m_data, current_scope, scope_path, sync_version] =
        m_thread_data.local();

    // Optional: Logging every stop might be spammy in parallel loops
    logger().trace(
        fmt::runtime(log_fmt_text), current_scope.to_string(), time_ms);

    assert(m_data.contains(current_scope));
    assert(m_data.at(current_scope).contains("time_ms"));
    assert(m_data.at(current_scope).contains("count"));
    m_data[current_scope]["time_ms"] =
        m_data[current_scope]["time_ms"].get<double>() + time_ms;
    m_data[current_scope]["count"] =
        m_data[current_scope]["count"].get<size_t>() + 1;

    if (std::this_thread::get_id() == m_main_thread_id) {
        size_t v = m_scope_version.load(std::memory_order_relaxed);
        m_scope_version.store(v + 1, std::memory_order_release);

        scope_path.pop_back();
        current_scope.pop_back();

        m_scope_version.store(v + 2, std::memory_order_release);
    } else {
        scope_path.pop_back();
        current_scope.pop_back();
    }
}

void Profiler::reset()
{
    m_thread_data.clear();
    m_scope_version = 0;
    m_main_thread_data = nullptr;
}

void Profiler::merge_json(nlohmann::json& target, const nlohmann::json& source)
{
    for (auto it = source.begin(); it != source.end(); ++it) {
        if (target.contains(it.key())) {
            if (it.value().is_number()) {
                // Sum numeric values (time_ms, count)
                if (target[it.key()].is_number()) {
                    if (it.value().is_number_integer()) {
                        target[it.key()] =
                            target[it.key()].get<int>() + it.value().get<int>();
                    } else if (it.value().is_number_float()) {
                        target[it.key()] = target[it.key()].get<double>()
                            + it.value().get<double>();
                    }
                }
            } else if (it.value().is_object()) {
                // Recursively merge objects (scopes)
                merge_json(target[it.key()], it.value());
            }
        } else {
            // New key, just copy
            target[it.key()] = it.value();
        }
    }
}

nlohmann::json Profiler::combine_data() const
{
    nlohmann::json combined = nlohmann::json::object();
    for (const auto& tld : m_thread_data) {
        merge_json(combined, tld.m_data);
    }
    return combined;
}

nlohmann::json Profiler::data() const { return combine_data(); }

void Profiler::print() const
{
    logger().info(
        "[{}] profiler: {}",
        fmt::format(fmt::fg(fmt::terminal_color::magenta), "timing"),
        combine_data().dump(2));
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
    nlohmann::json combined = combine_data();

    os << "Id,Parent,Name,Time (ms),Count\n";

    if (combined.empty()) {
        os << std::flush;
        return;
    }

    const nlohmann::json::json_pointer root;
    std::queue<std::pair<int, nlohmann::json::json_pointer>> queue;
    queue.push(std::make_pair(-1, root));
    int id = -1;

    while (!queue.empty()) {
        const auto [parent_id, ptr] = queue.front();
        queue.pop();

        assert(combined.contains(ptr));
        const auto& data = ptr == root ? combined : combined.at(ptr);
        if (ptr != root) {
            os << fmt::format(
                "{:d},{},{},{:.6g},{:d}\n", id,
                parent_id == -1 ? "" : std::to_string(parent_id), ptr.back(),
                data.at("time_ms").get<double>(),
                data.at("count").get<size_t>());
        }

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