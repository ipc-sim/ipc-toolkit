#include "distance_based_potential.hpp"

namespace ipc {

double DistanceBasedPotential::operator()(
    const Contact& contact, const VectorMax12d& x) const
{
    // w * m(x) * f(d(x))
    // NOTE: can save a multiplication by checking if !contact.is_mollified()
    return contact.weight * contact.mollifier(x)
        * distance_based_potential(contact.compute_distance(x), contact.dmin);
}

VectorMax12d DistanceBasedPotential::gradient(
    const Contact& contact, const VectorMax12d& x) const
{
    const double d = contact.compute_distance(x);                     // d(x)
    const VectorMax12d grad_d = contact.compute_distance_gradient(x); // ∇d(x)

    // f(d(x))
    const double f = distance_based_potential(d, contact.dmin);
    // f'(d(x))
    const double grad_f = distance_based_potential_gradient(d, contact.dmin);

    if (!contact.is_mollified()) {
        // ∇[f(d(x))] = f'(d(x)) * ∇d(x)
        return (contact.weight * grad_f) * grad_d;
    }

    const double m = contact.mollifier(x);                     // m(x)
    const VectorMax12d grad_m = contact.mollifier_gradient(x); // ∇m(x)

    // ∇[m(x) * f(d(x))] = f(d(x)) * ∇m(x) + m(x) * ∇ f(d(x))
    return (contact.weight * f) * grad_m
        + (contact.weight * m * grad_f) * grad_d;
}

MatrixMax12d DistanceBasedPotential::hessian(
    const Contact& contact,
    const VectorMax12d& x,
    const bool project_hessian_to_psd) const
{
    const double d = contact.compute_distance(x);                     // d(x)
    const VectorMax12d grad_d = contact.compute_distance_gradient(x); // ∇d(x)
    const MatrixMax12d hess_d = contact.compute_distance_hessian(x); // ∇²d(x)

    // f'(d(x))
    const double grad_f = distance_based_potential_gradient(d, contact.dmin);
    // f"(d(x))
    const double hess_f = distance_based_potential_hessian(d, contact.dmin);

    MatrixMax12d hess;
    if (!contact.is_mollified()) {
        // ∇²[f(d(x))] = ∇(f'(d(x)) * ∇d(x))
        //             = f"(d(x)) * ∇d(x) * ∇d(x)ᵀ + f'(d(x)) * ∇²d(x)
        hess = (contact.weight * hess_f) * grad_d * grad_d.transpose()
            + (contact.weight * grad_f) * hess_d;
    } else {
        const double f = distance_based_potential(d, contact.dmin); // f(d(x))

        const double m = contact.mollifier(x);                     // m(x)
        const VectorMax12d grad_m = contact.mollifier_gradient(x); // ∇ m(x)
        const MatrixMax12d hess_m = contact.mollifier_hessian(x);  // ∇² m(x)

        const double weighted_m = contact.weight * m;

        // ∇f(d(x)) * ∇m(x)ᵀ
        const MatrixMax12d grad_f_grad_m =
            (contact.weight * grad_f) * grad_d * grad_m.transpose();

        // ∇²[m(x) * b(d(x))] = ∇[∇m(x) * b(d(x)) + m(x) * ∇b(d(x))]
        //                    = ∇²m(x) * b(d(x)) + ∇b(d(x)) * ∇m(x)ᵀ
        //                      + ∇m(x) * ∇b(d(x))ᵀ + m(x) * ∇²b(d(x))
        hess = (contact.weight * f) * hess_m + grad_f_grad_m
            + grad_f_grad_m.transpose()
            + (weighted_m * hess_f) * grad_d * grad_d.transpose()
            + (weighted_m * grad_f) * hess_d;
    }

    // Need to project entire hessian because w can be negative
    return project_hessian_to_psd ? project_to_psd(hess) : hess;
}

} // namespace ipc