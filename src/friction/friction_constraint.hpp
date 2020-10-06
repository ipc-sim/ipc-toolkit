#pragma once

#include <ipc/collision_constraint.hpp>
#include <ipc/utils/eigen_ext.hpp>

namespace ipc {

struct FrictionConstraint {
    /// @brief Barycentric coordinates of the closest point(s)
    Eigen::VectorX2d closest_point;

    /// @brief Tangent basis of the contact (max size 3×2)
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor, 3, 2>
        tangent_basis;

    /// @brief Contact force magnitude
    double normal_force_magnitude;

    /// @brief Coefficient of friction
    double mu;
};

struct VertexVertexFrictionConstraint : VertexVertexConstraint,
                                        FrictionConstraint {
    VertexVertexFrictionConstraint(long vertex0_index, long vertex1_index);
    VertexVertexFrictionConstraint(const VertexVertexConstraint& constraint);
};

struct EdgeVertexFrictionConstraint : EdgeVertexConstraint, FrictionConstraint {
    EdgeVertexFrictionConstraint(long edge_index, long vertex_index);
    EdgeVertexFrictionConstraint(const EdgeVertexConstraint& constraint);
};

struct EdgeEdgeFrictionConstraint : EdgeEdgeConstraint, FrictionConstraint {
    EdgeEdgeFrictionConstraint(long edge0_index, long edge1_index);
    EdgeEdgeFrictionConstraint(const EdgeEdgeConstraint& constraint);
};

struct FaceVertexFrictionConstraint : FaceVertexConstraint, FrictionConstraint {
    FaceVertexFrictionConstraint(long face_index, long vertex_index);
    FaceVertexFrictionConstraint(const FaceVertexConstraint& constraint);
};

struct FrictionConstraints {
    std::vector<VertexVertexFrictionConstraint> vv_constraints;
    std::vector<EdgeVertexFrictionConstraint> ev_constraints;
    std::vector<EdgeEdgeFrictionConstraint> ee_constraints;
    std::vector<FaceVertexFrictionConstraint> fv_constraints;

    size_t size() const;

    size_t num_constraints() const;

    void clear();

    FrictionConstraint& operator[](size_t idx);
    const FrictionConstraint& operator[](size_t idx) const;
};

} // namespace ipc
