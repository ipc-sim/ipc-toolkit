#pragma once

#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

/// @brief A spatial hash for fast nearest neighbor queries.
/// @note Based on the implementation from Matthias Müller (https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/11-hashing.html)
class PointPointSpatialHash {
public:
    /// @brief Destroy the Spatial Hash object
    virtual ~PointPointSpatialHash() = default;

    /// @brief Build the spatial hash from a list of points.
    /// @param points Matrix of points (one per row).
    void build(const Eigen::MatrixXd& points);

    /// @brief Query the spatial hash for points within a maximum distance from a point.
    /// @param point The point to search around.
    /// @param max_dist The maximum distance to search.
    /// @return A list of point indices within the maximum distance from the point.
    const std::vector<int>&
    query(const VectorMax3d& point, const double max_dist) const;

protected:
    /// @brief Convert a coordinate to an integer coordinate.
    /// @param coords The coordinate to convert.
    /// @return The integer coordinate.
    virtual VectorMax3i int_coords(const VectorMax3d& coords) const = 0;

    /// @brief Hash a point to a cell index.
    /// @param coords Integer coordinates of the point.
    /// @return The cell index the point hashes to.
    virtual int hash(const VectorMax3i& coords) const = 0;

    /// @brief Hash a point to a cell index.
    /// @param point The point to hash.
    /// @return The cell index the point hashes to.
    int hash(const VectorMax3d& point) const;

    /// One over the spacing between cells in the spatial hash.
    double inv_spacing;
    /// Size of the spatial hash table.
    int table_size;
    /// Start index of each cell in the spatial hash.
    std::vector<int> cell_starts;
    /// List of point indices in the spatial hash.
    std::vector<int> cell_entries;

private:
    /// List of point indices in the spatial hash.
    mutable std::vector<int> query_ids;
};

class FinitePointPointSpatialHash : public PointPointSpatialHash {
public:
    FinitePointPointSpatialHash() = default;

    /// @brief Construct a new Spatial Hash object
    /// @param spacing The spacing between cells in the spatial hash.
    /// @param max_num_objects Maximum number of points that can be stored in the spatial hash.
    FinitePointPointSpatialHash(
        const VectorMax3d& domain_min,
        const VectorMax3d& domain_max,
        const double spacing,
        const int max_num_objects);

protected:
    /// @brief Convert a coordinate to an integer coordinate.
    /// @param coords The coordinate to convert.
    /// @return The integer coordinate.
    VectorMax3i int_coords(const VectorMax3d& coords) const override;

    /// @brief Hash a point to a cell index.
    /// @param coords Integer coordinates of the point.
    /// @return The cell index the point hashes to.
    int hash(const VectorMax3i& coords) const override;

    VectorMax3d domain_min;  ///< Minimum point in the spatial hash domain.
    VectorMax3i domain_size; ///< Number of cells in the x/y/z-directions.
};

class InfinitePointPointSpatialHash : public PointPointSpatialHash {
public:
    InfinitePointPointSpatialHash() = default;

    /// @brief Construct a new Spatial Hash object
    /// @param spacing The spacing between cells in the spatial hash.
    /// @param max_num_objects Maximum number of points that can be stored in the spatial hash.
    InfinitePointPointSpatialHash(
        const double spacing, const int max_num_objects);

protected:
    /// @brief Convert a coordinate to an integer coordinate.
    /// @param coords The coordinate to convert.
    /// @return The integer coordinate.
    VectorMax3i int_coords(const VectorMax3d& coords) const override;

    /// @brief Hash a point to a cell index.
    /// @param coords Integer coordinates of the point.
    /// @return The cell index the point hashes to.
    int hash(const VectorMax3i& coords) const override;
};

} // namespace ipc