#pragma once

#include <ipc/candidates/collision_stencil.hpp>
#include <ipc/utils/eigen_ext.hpp>
#include <ipc/collision_mesh.hpp>

#include <Eigen/Core>

#include <array>

namespace ipc {

class NormalCollision : virtual public CollisionStencil {
public:
    NormalCollision() = default;

    NormalCollision(
        const double weight,
        const Eigen::SparseVector<double>& weight_gradient);

    virtual ~NormalCollision() = default;

    // -- Distance mollifier ---------------------------------------------------

    /// @brief Does the distance potentially have to be mollified?
    virtual bool is_mollified() const { return false; }

    /// @brief Compute the mollifier threshold for the distance.
    /// @param rest_positions The stencil's rest vertex positions.
    /// @return The mollifier threshold.
    virtual double
    mollifier_threshold(Eigen::ConstRef<VectorMax12d> rest_positions) const
    {
        return std::numeric_limits<double>::quiet_NaN(); // No mollifier
    }

    /// @brief Compute the mollifier for the distance.
    /// @param positions The stencil's vertex positions.
    /// @return The mollifier value.
    virtual double mollifier(Eigen::ConstRef<VectorMax12d> positions) const;

    /// @brief Compute the mollifier for the distance.
    /// @param positions The stencil's vertex positions.
    /// @param eps_x The mollifier's threshold.
    /// @return The mollifier value.
    virtual double
    mollifier(Eigen::ConstRef<VectorMax12d> positions, double eps_x) const;

    /// @brief Compute the gradient of the mollifier for the distance wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @return The mollifier gradient.
    virtual VectorMax12d
    mollifier_gradient(Eigen::ConstRef<VectorMax12d> positions) const;

    /// @brief Compute the gradient of the mollifier for the distance wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @param eps_x The mollifier's threshold.
    /// @return The mollifier gradient.
    virtual VectorMax12d mollifier_gradient(
        Eigen::ConstRef<VectorMax12d> positions, double eps_x) const;

    /// @brief Compute the Hessian of the mollifier for the distance wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @return The mollifier Hessian.
    virtual MatrixMax12d
    mollifier_hessian(Eigen::ConstRef<VectorMax12d> positions) const;

    /// @brief Compute the Hessian of the mollifier for the distance wrt the positions.
    /// @param positions The stencil's vertex positions.
    /// @param eps_x The mollifier's threshold.
    /// @return The mollifier Hessian.
    virtual MatrixMax12d mollifier_hessian(
        Eigen::ConstRef<VectorMax12d> positions, double eps_x) const;

    /// @brief Compute the gradient of the mollifier for the distance w.r.t. rest positions.
    /// @param rest_positions The stencil's rest vertex positions.
    /// @param positions The stencil's vertex positions.
    /// @return The mollifier gradient w.r.t. rest positions.
    virtual Vector12d mollifier_gradient_wrt_x(
        Eigen::ConstRef<VectorMax12d> rest_positions,
        Eigen::ConstRef<VectorMax12d> positions) const;

    /// @brief Compute the jacobian of the distance mollifier's gradient w.r.t. rest positions.
    /// @param rest_positions The stencil's rest vertex positions.
    /// @param positions The stencil's vertex positions.
    /// @return The jacobian of the mollifier's gradient w.r.t. rest positions.
    virtual Matrix12d mollifier_gradient_jacobian_wrt_x(
        Eigen::ConstRef<VectorMax12d> rest_positions,
        Eigen::ConstRef<VectorMax12d> positions) const;

    // -------------------------------------------------------------------------

    /// @brief The minimum separation distance.
    double dmin = 0;

    /// @brief The term's weight (e.g., collision area)
    double weight = 1;

    /// @brief The gradient of the term's weight wrt the rest positions.
    Eigen::SparseVector<double> weight_gradient;

    /// @brief Material ID for the first object in the collision.
    /// Set to NO_MATERIAL_ID if material-specific friction is not needed.
    int material_id1 = NO_MATERIAL_ID;

    /// @brief Material ID for the second object in the collision.
    /// Set to NO_MATERIAL_ID if material-specific friction is not needed.
    int material_id2 = NO_MATERIAL_ID;

    /// @brief Check if material IDs are being used for this collision
    /// @return true if material IDs are being used, false otherwise
    bool has_material_ids() const { return material_id1 != NO_MATERIAL_ID && material_id2 != NO_MATERIAL_ID; }
    
    /// @brief Set the material IDs for this collision
    /// @param mat_id1 Material ID for the first object
    /// @param mat_id2 Material ID for the second object
    void set_material_ids(int mat_id1, int mat_id2) 
    {
        material_id1 = mat_id1;
        material_id2 = mat_id2;
    }

    /// @brief Get the material IDs involved in this collision
    /// @return Pair of material IDs (first, second)
    std::pair<int, int> get_material_ids() const
    {
        return {material_id1, material_id2};
    }
};

} // namespace ipc
