#pragma once

namespace ipc::autogen {

// hess is (144×1) flattened in column-major order
void edge_edge_closest_point_hessian_a(
    double ea0_x,
    double ea0_y,
    double ea0_z,
    double ea1_x,
    double ea1_y,
    double ea1_z,
    double eb0_x,
    double eb0_y,
    double eb0_z,
    double eb1_x,
    double eb1_y,
    double eb1_z,
    double hess[144]);

// hess is (144×1) flattened in column-major order
void edge_edge_closest_point_hessian_b(
    double ea0_x,
    double ea0_y,
    double ea0_z,
    double ea1_x,
    double ea1_y,
    double ea1_z,
    double eb0_x,
    double eb0_y,
    double eb0_z,
    double eb1_x,
    double eb1_y,
    double eb1_z,
    double hess[144]);

void face_normal_squared_norm_gradient(
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double grad[9]);

// hess is (81×1) flattened in column-major order
void face_normal_squared_norm_hessian(
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double hess[81]);

void face_term_aux_gradient(
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double p1_x,
    double p1_y,
    double p1_z,
    double p2_x,
    double p2_y,
    double p2_z,
    double grad[15]);

// hess is (225×1) flattened in column-major order
void face_term_aux_hessian(
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double p1_x,
    double p1_y,
    double p1_z,
    double p2_x,
    double p2_y,
    double p2_z,
    double hess[225]);

// hess is (144×1) flattened in column-major order
void triangle_closest_point_hessian_0(
    double p_x,
    double p_y,
    double p_z,
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double hess[144]);
// hess is (144×1) flattened in column-major order
void triangle_closest_point_hessian_1(
    double p_x,
    double p_y,
    double p_z,
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double hess[144]);

void face_term_aux_fast_gradient(
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double p_x,
    double p_y,
    double p_z,
    double d,
    double grad[13]);

// hess is (169×1) flattened in column-major order
void face_term_aux_fast_hessian(
    double t0_x,
    double t0_y,
    double t0_z,
    double t1_x,
    double t1_y,
    double t1_z,
    double t2_x,
    double t2_y,
    double t2_z,
    double p_x,
    double p_y,
    double p_z,
    double d,
    double hess[169]);

} // namespace ipc::autogen
