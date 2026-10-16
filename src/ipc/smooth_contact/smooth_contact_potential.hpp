#pragma once

// Deprecated: include <ipc/gcp/gcp_potential.hpp> instead. We warn at include
// time because the old header path has no symbol for [[deprecated]] to attach
// to. The compiler reports this header's path along with the message.
#if defined(_MSC_VER)
#pragma message(                                                               \
    "warning: ipc/smooth_contact/smooth_contact_potential.hpp is deprecated. " \
    "Include <ipc/gcp/gcp_potential.hpp> instead.")
#else
#pragma GCC warning                                                            \
    "Deprecated header. Include <ipc/gcp/gcp_potential.hpp> instead."
#endif

#include <ipc/gcp/gcp_potential.hpp>
