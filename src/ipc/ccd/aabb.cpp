#include "aabb.hpp"

namespace ipc {

// Discrete collision detection

bool point_edge_aabb_cd(
    const VectorMax3d& p,
    const VectorMax3d& e0,
    const VectorMax3d& e1,
    double dist)
{
    const ArrayMax3d max_e = e0.array().max(e1.array());
    const ArrayMax3d min_e = e0.array().min(e1.array());
    return (p.array() <= max_e + dist).all()
        && (min_e <= p.array() + dist).all();
}

bool edge_edge_aabb_cd(
    const VectorMax3d& ea0,
    const VectorMax3d& ea1,
    const VectorMax3d& eb0,
    const VectorMax3d& eb1,
    double dist)
{
    const ArrayMax3d max_a = ea0.array().max(ea1.array());
    const ArrayMax3d min_a = ea0.array().min(ea1.array());
    const ArrayMax3d max_b = eb0.array().max(eb1.array());
    const ArrayMax3d min_b = eb0.array().min(eb1.array());
    return (min_a <= max_b + dist).all() && (min_b <= max_a + dist).all();
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
    return (p.array() <= max_tri + dist).all()
        && (min_tri <= p.array() + dist).all();
}

bool edge_triangle_aabb_cd(
    const Eigen::Vector3d& e0,
    const Eigen::Vector3d& e1,
    const Eigen::Vector3d& t0,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2,
    double dist)
{
    const Eigen::Array3d max_e = e0.array().max(e1.array());
    const Eigen::Array3d min_e = e0.array().min(e1.array());
    const Eigen::Array3d max_tri = t0.array().max(t1.array()).max(t2.array());
    const Eigen::Array3d min_tri = t0.array().min(t1.array()).min(t2.array());
    return (min_e <= max_tri + dist).all() && (min_tri <= max_e + dist).all();
}

// Continous collision detection

bool point_edge_aabb_ccd(
    const VectorMax3d& p_t0,
    const VectorMax3d& e0_t0,
    const VectorMax3d& e1_t0,
    const VectorMax3d& p_t1,
    const VectorMax3d& e0_t1,
    const VectorMax3d& e1_t1,
    double dist)
{
    const ArrayMax3d max_p = p_t0.array().max(p_t1.array());
    const ArrayMax3d min_p = p_t0.array().min(p_t1.array());
    const ArrayMax3d max_e =
        e0_t0.array().max(e1_t0.array()).max(e0_t1.array()).max(e1_t1.array());
    const ArrayMax3d min_e =
        e0_t0.array().min(e1_t0.array()).min(e0_t1.array()).min(e1_t1.array());
    return (min_p <= max_e + dist).all() && (min_e <= max_p + dist).all();
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
    return (min_p <= max_tri + dist).all() && (min_tri <= max_p + dist).all();
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
    return (min_a <= max_b + dist).all() && (min_b <= max_a + dist).all();
}

} // namespace ipc
