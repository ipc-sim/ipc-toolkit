
#include <npe.h>

#include <ipc/barrier/barrier.hpp>

#include <common.hpp>

#define NPE_FUNCTION(name) npe_function(name)

const char* ds_barrier = R"ipc_Qu8mg5v7(
Function that grows to infinity as d approaches 0 from the right.

Parameters
----------
d : The distance
dhat : Activation distance of the barrier

Returns
-------
The value of the barrier function at d.

See also
--------

Notes
-----
.. math:: b(d) = -(d-\hat{d})^2\ln(\frac{d}{\hat{d}})

Examples
--------

)ipc_Qu8mg5v7";

// clang-format off
npe_function(barrier)
npe_doc(ds_barrier)
npe_arg(d, double)
npe_arg(dhat, double)
npe_begin_code()

    return ipc::barrier(d, dhat);

npe_end_code()
; // clang-format on
