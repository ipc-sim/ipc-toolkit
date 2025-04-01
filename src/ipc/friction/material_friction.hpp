#pragma once

#include <ipc/collisions/tangential/tangential_collision.hpp>
#include <Eigen/Core>
#include <memory>

namespace ipc {

/// @brief A class to handle material-specific friction coefficients
class MaterialFriction {
public:
    /// @brief Constructor
    MaterialFriction() = default;
    
    /// @brief Constructor with a friction table
    /// @param friction_table Table of friction coefficients indexed by material IDs
    MaterialFriction(const Eigen::MatrixXd& friction_table)
        : m_friction_table(std::make_shared<Eigen::MatrixXd>(friction_table))
    {}
    
    /// @brief Set the friction coefficient table
    /// @param friction_table Table of friction coefficients indexed by material IDs
    void set_friction_table(const Eigen::MatrixXd& friction_table)
    {
        m_friction_table = std::make_shared<Eigen::MatrixXd>(friction_table);
    }
    
    /// @brief Get the friction coefficient table
    /// @return Shared pointer to the friction table
    std::shared_ptr<const Eigen::MatrixXd> get_friction_table() const
    {
        return m_friction_table;
    }
    
    /// @brief Get the friction coefficient for a collision
    /// @param collision The tangential collision
    /// @param default_mu The default friction coefficient to use if no material IDs are set
    /// @return The friction coefficient based on material IDs, or default_mu if no match
    double get_coefficient(const TangentialCollision& collision, double default_mu) const
    {
        // If material IDs are not set or there's no table, use default
        if (!collision.has_material_ids() || !m_friction_table || m_friction_table->size() == 0) {
            return default_mu;
        }
        
        int id1 = collision.material_id1;
        int id2 = collision.material_id2;
        
        // Check if material IDs are valid for the friction table
        if (id1 >= 0 && id1 < m_friction_table->rows() && 
            id2 >= 0 && id2 < m_friction_table->cols()) {
            return (*m_friction_table)(id1, id2);
        }
        
        // Default to collision's mu if material IDs are out of range
        return default_mu;
    }
    
private:
    /// @brief Table of friction coefficients indexed by material IDs
    std::shared_ptr<Eigen::MatrixXd> m_friction_table;
};

} // namespace ipc
