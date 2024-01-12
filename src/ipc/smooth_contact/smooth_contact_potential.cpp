#include "smooth_contact_potential.hpp"

namespace ipc {

// template class SmoothContactPotential<VirtualCollisions<4>>;
template class SmoothContactPotential<SmoothCollisions<2>>;
template class SmoothContactPotential<SmoothCollisions<3>>;

} // namespace ipc