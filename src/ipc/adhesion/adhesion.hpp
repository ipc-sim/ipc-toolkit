#pragma once

namespace ipc {

// Fang and Li et al. [2023]:

double normal_adhesion_potential(
    const double d, const double dhat_p, const double dhat_a, const double a2);

double normal_adhesion_potential_gradient(
    const double d, const double dhat_p, const double dhat_a, const double a2);

double normal_adhesion_potential_hessian(
    const double d, const double dhat_p, const double dhat_a, const double a2);

} // namespace ipc
