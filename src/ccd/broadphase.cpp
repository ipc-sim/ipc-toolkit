#include <ipc/ccd/broadphase.hpp>

namespace ipc {

// Discrete collision detection

bool point_edge_aabb_cd(
    const Eigen::Vector2d& p,
    const Eigen::Vector2d& e0,
    const Eigen::Vector2d& e1,
    double dist)
{
    const Eigen::Array2d max_e = e0.array().max(e1.array());
    const Eigen::Array2d min_e = e0.array().min(e1.array());
    if ((p.array() > max_e + dist).any() || (min_e > p.array() + dist).any()) {
        return false;
    } else {
        return true;
    }
}

bool point_triangle_aabb_cd(
    const Eigen::Vector3d& p,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2,
    double dist)
{
    const Eigen::Array3d max_tri = t0.array().max(t1.array()).max(t2.array());
    const Eigen::Array3d min_tri = t0.array().min(t1.array()).min(t2.array());
    if ((p.array() > max_tri + dist).any()
        || (min_tri > p.array() + dist).any()) {
        return false;
    } else {
        return true;
    }
}

bool edge_edge_aabb_cd(
    const Eigen::Vector3d& ea0,
    const Eigen::Vector3d& ea1,
    const Eigen::Vector3d& eb0,
    const Eigen::Vector3d& eb1,
    double dist)
{
    const Eigen::Array3d max_a = ea0.array().max(ea1.array());
    const Eigen::Array3d min_a = ea0.array().min(ea1.array());
    const Eigen::Array3d max_b = eb0.array().max(eb1.array());
    const Eigen::Array3d min_b = eb0.array().min(eb1.array());
    if ((min_a > max_b + dist).any() || (min_b > max_a + dist).any()) {
        return false;
    } else {
        return true;
    }
}

// Continous collision detection

bool point_edge_aabb_ccd(
    const Eigen::VectorX3d& p_t0,
    const Eigen::VectorX3d& e0_t0,
    const Eigen::VectorX3d& e1_t0,
    const Eigen::VectorX3d& p_t1,
    const Eigen::VectorX3d& e0_t1,
    const Eigen::VectorX3d& e1_t1,
    double dist)
{
    const Eigen::ArrayMax3d max_p = p_t0.array().max(p_t1.array());
    const Eigen::ArrayMax3d min_p = p_t0.array().min(p_t1.array());
    const Eigen::ArrayMax3d max_e =
        e0_t0.array().max(e1_t0.array()).max(e0_t1.array()).max(e1_t1.array());
    const Eigen::ArrayMax3d min_e =
        e0_t0.array().min(e1_t0.array()).min(e0_t1.array()).min(e1_t1.array());
    if ((min_p > max_e + dist).any() || (min_e > max_p + dist).any()) {
        return false;
    } else {
        return true;
    }
}

bool point_triangle_aabb_ccd(
    const Eigen::Vector3d& p_t0,
    const Eigen::Vector3d& t0_t0,
    const Eigen::Vector3d& t1_t0,
    const Eigen::Vector3d& t2_t0,
    const Eigen::Vector3d& p_t1,
    const Eigen::Vector3d& t0_t1,
    const Eigen::Vector3d& t1_t1,
    const Eigen::Vector3d& t2_t1,
    double dist)
{
    const Eigen::Array3d max_p = p_t0.array().max((p_t1).array());
    const Eigen::Array3d min_p = p_t0.array().min((p_t1).array());
    const Eigen::Array3d max_tri = t0_t0.array()
                                       .max(t1_t0.array())
                                       .max(t2_t0.array())
                                       .max(t0_t1.array())
                                       .max(t1_t1.array())
                                       .max(t2_t1.array());
    const Eigen::Array3d min_tri = t0_t0.array()
                                       .min(t1_t0.array())
                                       .min(t2_t0.array())
                                       .min(t0_t1.array())
                                       .min(t1_t1.array())
                                       .min(t2_t1.array());
    if ((min_p > max_tri + dist).any() || (min_tri > max_p + dist).any()) {
        return false;
    } else {
        return true;
    }
}

bool edge_edge_aabb_ccd(
    const Eigen::Vector3d& ea0_t0,
    const Eigen::Vector3d& ea1_t0,
    const Eigen::Vector3d& eb0_t0,
    const Eigen::Vector3d& eb1_t0,
    const Eigen::Vector3d& ea0_t1,
    const Eigen::Vector3d& ea1_t1,
    const Eigen::Vector3d& eb0_t1,
    const Eigen::Vector3d& eb1_t1,
    double dist)
{
    const Eigen::Array3d max_a = ea0_t0.array()
                                     .max(ea1_t0.array())
                                     .max(ea0_t1.array())
                                     .max(ea1_t1.array());
    const Eigen::Array3d min_a = ea0_t0.array()
                                     .min(ea1_t0.array())
                                     .min(ea0_t1.array())
                                     .min(ea1_t1.array());
    const Eigen::Array3d max_b = eb0_t0.array()
                                     .max(eb1_t0.array())
                                     .max(eb0_t1.array())
                                     .max(eb1_t1.array());
    const Eigen::Array3d min_b = eb0_t0.array()
                                     .min(eb1_t0.array())
                                     .min(eb0_t1.array())
                                     .min(eb1_t1.array());
    if ((min_a > max_b + dist).any() || (min_b > max_a + dist).any()) {
        return false;
    } else {
        return true;
    }
}

} // namespace ipc
