// clang-format off

#include <npe.h>

#include <ipc/distance/point_point.hpp>

#include <common.hpp>

const char* ds_point_point_distance = R"ipc_Qu8mg5v7(
Compute the distance between two points.

Parameters
----------
p0 : The first point
p1 : The second point

Returns
-------
The distance between p0 and p1

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_point_distance)
npe_doc(ds_point_point_distance)
npe_arg(p0, dense_double)
npe_arg(p1, npe_matches(p0))
npe_begin_code()

    assert_2D_or_3D_vector(p0, "p0");
    assert_2D_or_3D_vector(p1, "p1");

    return ipc::point_point_distance(p0, p1);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_point_distance_gradient = R"ipc_Qu8mg5v7(
Compute the gradient of the distance between two points.

Parameters
----------
p0 : The first point
p1 : The second point

Returns
-------
The gradient of the distance wrt p0 and p1.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_point_distance_gradient)
npe_doc(ds_point_point_distance_gradient)
npe_arg(p0, dense_double)
npe_arg(p1, npe_matches(p0))
npe_begin_code()

    assert_2D_or_3D_vector(p0, "p0");
    assert_2D_or_3D_vector(p1, "p1");

    Eigen::VectorX12<npe_Scalar_p0> grad;
    ipc::point_point_distance_gradient(p0, p1, grad);
    return npe::move(grad);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_point_distance_hessian = R"ipc_Qu8mg5v7(
Compute the hessian of the distance between a point and point.

Parameters
----------
p : The point
p0 : The first vertex of the point
p1 : The second vertex of the point
dtype : (Optional) The point point distance type to compute

Returns
-------
The hessian of the distance wrt p, p0, and p1.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_point_distance_hessian)
npe_doc(ds_point_point_distance_hessian)
npe_arg(p0, dense_double)
npe_arg(p1, npe_matches(p0))
npe_begin_code()

    assert_2D_or_3D_vector(p0, "p0");
    assert_2D_or_3D_vector(p1, "p1");

    Eigen::MatrixXX12<npe_Scalar_p0> hess;
    ipc::point_point_distance_hessian(p0, p1, hess);
    return npe::move(hess);

npe_end_code()
