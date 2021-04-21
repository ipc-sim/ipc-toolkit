// clang-format off

#include <npe.h>

#include <ipc/distance/edge_edge_mollifier.hpp>

#include <common.hpp>

const char* ds_edge_edge_cross_squarednorm = R"ipc_Qu8mg5v7(
Compute the squared norm of the edge-edge cross product.

Parameters
----------
ea0 : The first vertex of the first edge
ea1 : The second vertex of the first edge
eb0 : The first vertex of the second edge
eb1 : The second vertex of the second edge

Returns
-------
The squared norm of the edge-edge cross product.

See also
--------

Notes
-----

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_cross_squarednorm)
npe_doc(ds_edge_edge_cross_squarednorm)
npe_arg(ea0, dense_double)
npe_arg(ea1, npe_matches(ea0))
npe_arg(eb0, npe_matches(ea0))
npe_arg(eb1, npe_matches(ea0))
npe_begin_code()

    assert_3D_vector(ea0, "ea0");
    assert_3D_vector(ea1, "ea1");
    assert_3D_vector(eb0, "eb0");
    assert_3D_vector(eb1, "eb1");

    Eigen::Vector3<npe_Scalar_ea0> ea0_copy, ea1_copy, eb0_copy, eb1_copy;
    copy_vector(ea0, ea0_copy);
    copy_vector(ea1, ea1_copy);
    copy_vector(eb0, eb0_copy);
    copy_vector(eb1, eb1_copy);

    return ipc::edge_edge_cross_squarednorm(ea0_copy, ea1_copy, eb0_copy, eb1_copy);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_cross_squarednorm_gradient = R"ipc_Qu8mg5v7(
Compute the gradient of the squared norm of the edge cross product.

Parameters
----------
ea0 : The first vertex of the first edge
ea1 : The second vertex of the first edge
eb0 : The first vertex of the second edge
eb1 : The second vertex of the second edge

Returns
-------
The gradient of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_cross_squarednorm_gradient)
npe_doc(ds_edge_edge_cross_squarednorm_gradient)
npe_arg(ea0, dense_double)
npe_arg(ea1, npe_matches(ea0))
npe_arg(eb0, npe_matches(ea0))
npe_arg(eb1, npe_matches(ea0))
npe_begin_code()

    assert_3D_vector(ea0, "ea0");
    assert_3D_vector(ea1, "ea1");
    assert_3D_vector(eb0, "eb0");
    assert_3D_vector(eb1, "eb1");

    Eigen::Vector3<npe_Scalar_ea0> ea0_copy, ea1_copy, eb0_copy, eb1_copy;
    copy_vector(ea0, ea0_copy);
    copy_vector(ea1, ea1_copy);
    copy_vector(eb0, eb0_copy);
    copy_vector(eb1, eb1_copy);

    Eigen::Vector<npe_Scalar_ea0, 12> grad;
    ipc::edge_edge_cross_squarednorm_gradient(
        ea0_copy, ea1_copy, eb0_copy, eb1_copy, grad);
    return npe::move(grad);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_cross_squarednorm_hessian = R"ipc_Qu8mg5v7(
Compute the hessian of the squared norm of the edge cross product.

Parameters
----------
ea0 : The first vertex of the first edge
ea1 : The second vertex of the first edge
eb0 : The first vertex of the second edge
eb1 : The second vertex of the second edge

Returns
-------
The hessian of the squared norm of the edge cross product wrt ea0, ea1, eb0, and eb1.

See also
--------

Notes
-----

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_cross_squarednorm_hessian)
npe_doc(ds_edge_edge_cross_squarednorm_hessian)
npe_arg(ea0, dense_double)
npe_arg(ea1, npe_matches(ea0))
npe_arg(eb0, npe_matches(ea0))
npe_arg(eb1, npe_matches(ea0))
npe_begin_code()

    assert_3D_vector(ea0, "ea0");
    assert_3D_vector(ea1, "ea1");
    assert_3D_vector(eb0, "eb0");
    assert_3D_vector(eb1, "eb1");

    Eigen::Vector3<npe_Scalar_ea0> ea0_copy, ea1_copy, eb0_copy, eb1_copy;
    copy_vector(ea0, ea0_copy);
    copy_vector(ea1, ea1_copy);
    copy_vector(eb0, eb0_copy);
    copy_vector(eb1, eb1_copy);

    Eigen::Matrix<npe_Scalar_ea0, 12, 12> hess;
    ipc::edge_edge_cross_squarednorm_hessian(
        ea0_copy, ea1_copy, eb0_copy, eb1_copy, hess);
    return npe::move(hess);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_mollifier = R"ipc_Qu8mg5v7(
Mollifier function for edge-edge distance.

Parameters
----------
x : Squared norm of the edge-edge cross product
eps_x : Mollifier activation threshold.

Returns
-------
The mollifier coefficient to premultiply the edge-edge distance.

See also
--------

Notes
-----

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_mollifier)
npe_doc(ds_edge_edge_mollifier)
npe_arg(x, double)
npe_arg(eps_x, double)
npe_begin_code()

    return ipc::edge_edge_mollifier(x, eps_x);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_mollifier_gradient = R"ipc_Qu8mg5v7(
The gradient of the mollifier function for edge-edge distance.

Parameters
----------
x : Squared norm of the edge-edge cross product
eps_x : Mollifier activation threshold.

Returns
-------
The gradient of the mollifier function for edge-edge distance wrt x.

See also
--------

Notes
-----

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_mollifier_gradient)
npe_doc(ds_edge_edge_mollifier_gradient)
npe_arg(x, double)
npe_arg(eps_x, double)
npe_begin_code()

    return ipc::edge_edge_mollifier_gradient(x, eps_x);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_mollifier_hessian = R"ipc_Qu8mg5v7(
The hessian of the mollifier function for edge-edge distance.

Parameters
----------
x : Squared norm of the edge-edge cross product
eps_x : Mollifier activation threshold.

Returns
-------
The hessian of the mollifier function for edge-edge distance wrt x.

See also
--------

Notes
-----

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_mollifier_hessian)
npe_doc(ds_edge_edge_mollifier_hessian)
npe_arg(x, double)
npe_arg(eps_x, double)
npe_begin_code()

    return ipc::edge_edge_mollifier_hessian(x, eps_x);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_mollifier_threshold = R"ipc_Qu8mg5v7(
Compute the threshold of the mollifier edge-edge distance.

This values is computed based on the edges at rest length.

Parameters
----------
x : Squared norm of the edge-edge cross product
eps_x : Mollifier activation threshold.

Returns
-------
Threshold for edge-edge mollification.

See also
--------

Notes
-----

Examples
--------

)ipc_Qu8mg5v7";
// clang-format on
NPE_FUNCTION(edge_edge_mollifier_threshold)
NPE_DOC(ds_edge_edge_mollifier_threshold)
NPE_ARG(ea0_rest, dense_double)
NPE_ARG(ea1_rest, NPE_MATCHES(ea0_rest))
NPE_ARG(eb0_rest, NPE_MATCHES(ea0_rest))
NPE_ARG(eb1_rest, NPE_MATCHES(ea0_rest))
NPE_BEGIN_CODE()

return ipc::edge_edge_mollifier_threshold(
    ea0_rest, ea1_rest, eb0_rest, eb1_rest);

NPE_END_CODE()
