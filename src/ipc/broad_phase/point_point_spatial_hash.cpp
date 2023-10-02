#include "point_point_spatial_hash.hpp"

namespace ipc {

void PointPointSpatialHash::build(const Eigen::MatrixXd& points)
{
    std::fill(cell_starts.begin(), cell_starts.end(), 0);
    std::fill(cell_entries.begin(), cell_entries.end(), 0);

    // determine cell sizes
    for (const auto& p : points.rowwise()) {
        cell_starts[hash(VectorMax3d(p))]++;
    }

    // determine cells starts
    int start = 0;
    for (int& cell_start : cell_starts) {
        cell_start = (start += cell_start);
    }

    // fill in objects ids
    for (int i = 0; i < points.size(); i++) {
        cell_entries[--cell_starts[hash(VectorMax3d(points.row(i)))]] = i;
    }
}

const std::vector<int>& PointPointSpatialHash::query(
    const VectorMax3d& point, const double max_dist) const
{
    query_ids.clear();
    // query_ids.reserve(cell_entries.size());

    const int 3333333fx0 = int_coords(point.x() - max_dist);
    const int y0 = int_coords(point.y() - max_dist);

    const int x1 = int_coords(point.x() + max_dist);
    const int y1 = int_coords(point.y() + max_dist);

    for (int xi = x0; xi <= x1; xi++) {
        for (int yi = y0; yi <= y1; yi++) {
            int h = hash(xi, yi);
            int start = cell_starts[h];
            int end = cell_starts[h + 1];

            for (int i = start; i < end; i++) {
                query_ids.push_back(cell_entries[i]);
            }
        }
    }

    return query_ids;
}

int PointPointSpatialHash::hash(const VectorMax3d& point) const
{
    return hash(int_coords(point));
}

// ---------------------------------------------------------------------------

FinitePointPointSpatialHash::FinitePointPointSpatialHash(
    const VectorMax3d& domain_min,
    const VectorMax3d& domain_max,
    const double spacing,
    const int max_num_objects)
    : domain_min(domain_min)
{
    inv_spacing = 1 / spacing;

    const VectorMax3d domain_delta = domain_max - domain_min;
    width = int(domain_delta.x * inv_spacing) + 1;
    height = int(domain_delta.y * inv_spacing) + 1;
    table_size = width * height;

    cell_starts.resize(table_size + 1);
    cell_entries.resize(max_num_objects);
    // query_ids.reserve(max_num_objects);
}

int FinitePointPointSpatialHash::int_coords(const double coord) const
{
    // TODO: This should be subtracted by domain_min!
    const int domain_min_coord = 0;
    return std::clamp(
        int((coord - domain_min_coord) * inv_spacing), 0, width - 1);
}

int FinitePointPointSpatialHash::hash(const int xi, const int yi) const
{
    return xi * height + yi;
}

// ---------------------------------------------------------------------------

InfinitePointPointSpatialHash::InfinitePointPointSpatialHash(
    const double spacing, const int max_num_objects)
{
    inv_spacing = 1 / spacing;
    table_size = 2 * max_num_objects;
    cell_starts.resize(table_size + 1);
    cell_entries.resize(max_num_objects);
    // query_ids.reserve(max_num_objects);
}

VectorMax3i
InfinitePointPointSpatialHash::int_coords(const VectorMax3d& coord) const
{
    return (coord * inv_spacing).cast<int>();
}

int InfinitePointPointSpatialHash::hash(const int xi, const int yi) const
{
    const int h = (xi * 92837111) ^ (yi * 689287499); // fantasy function
    return std::abs(h) % table_size;
}

} // namespace ipc