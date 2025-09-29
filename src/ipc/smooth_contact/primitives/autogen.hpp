#pragma once

namespace ipc {
namespace autogen {

    void edge_normal_term_gradient(
        double d_x,
        double d_y,
        double d_z,
        double e0_x,
        double e0_y,
        double e0_z,
        double e1_x,
        double e1_y,
        double e1_z,
        double f0_x,
        double f0_y,
        double f0_z,
        double f1_x,
        double f1_y,
        double f1_z,
        double grad[15]);
    // hess is (225×1) flattened in column-major order
    void edge_normal_term_hessian(
        double d_x,
        double d_y,
        double d_z,
        double e0_x,
        double e0_y,
        double e0_z,
        double e1_x,
        double e1_y,
        double e1_z,
        double f0_x,
        double f0_y,
        double f0_z,
        double f1_x,
        double f1_y,
        double f1_z,
        double hess[225]);
} // namespace autogen
} // namespace ipc
