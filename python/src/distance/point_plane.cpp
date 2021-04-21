// clang-format off

#include <npe.h>

#include <ipc/distance/point_plane.hpp>

#include <common.hpp>

const char* ds_point_plane_distance = R"ipc_Qu8mg5v7(
Compute the distance between a point and a plane.

Parameters
----------
p  : The point.
t0 : The first vertex of the triangle.
t1 : The second vertex of the triangle.
t2 : The third vertex of the triangle.

Returns
-------
The distance between the point and plane.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_plane_distance)
npe_doc(ds_point_plane_distance)
npe_arg(p, dense_double)
npe_arg(t0, npe_matches(p))
npe_arg(t1, npe_matches(p))
npe_arg(t2, npe_matches(p))
npe_begin_code()

    assert_3D_vector(p, "p");
    assert_3D_vector(t0, "t0");
    assert_3D_vector(t1, "t1");
    assert_3D_vector(t2, "t2");

    Eigen::Vector3<npe_Scalar_p> p_copy, t0_copy, t1_copy, t2_copy;
    copy_vector(p, p_copy);
    copy_vector(t0, t0_copy);
    copy_vector(t1, t1_copy);
    copy_vector(t2, t2_copy);

    return ipc::point_plane_distance(p_copy, t0_copy, t1_copy, t2_copy);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_plane_distance_gradient = R"ipc_Qu8mg5v7(
Compute the gradient of the distance between a point and a plane.

Parameters
----------
p  : The point
t0 : The first vertex of the triangle
t1 : The second vertex of the triangle
t2 : The third vertex of the triangle

Returns
-------
The gradient of the distance wrt p, t0, t1, and t2.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_plane_distance_gradient)
npe_doc(ds_point_plane_distance_gradient)
npe_arg(p, dense_double)
npe_arg(t0, npe_matches(p))
npe_arg(t1, npe_matches(p))
npe_arg(t2, npe_matches(p))
npe_begin_code()

    assert_3D_vector(p, "p");
    assert_3D_vector(t0, "t0");
    assert_3D_vector(t1, "t1");
    assert_3D_vector(t2, "t2");

    Eigen::Vector3<npe_Scalar_p> p_copy, t0_copy, t1_copy, t2_copy;
    copy_vector(p, p_copy);
    copy_vector(t0, t0_copy);
    copy_vector(t1, t1_copy);
    copy_vector(t2, t2_copy);

    Eigen::Vector<npe_Scalar_p, 12> grad;
    ipc::point_plane_distance_gradient(p_copy, t0_copy, t1_copy, t2_copy, grad);
    return npe::move(grad);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_plane_distance_hessian = R"ipc_Qu8mg5v7(
Compute the hessian of the distance between a point and a plane.

Parameters
----------
p  : The point
t0 : The first vertex of the triangle
t1 : The second vertex of the triangle
t2 : The third vertex of the triangle

Returns
-------
The hessian of the distance wrt p, t0, t1, and t2.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_plane_distance_hessian)
npe_doc(ds_point_plane_distance_hessian)
npe_arg(p, dense_double)
npe_arg(t0, npe_matches(p))
npe_arg(t1, npe_matches(p))
npe_arg(t2, npe_matches(p))
npe_begin_code()

    assert_3D_vector(p, "p");
    assert_3D_vector(t0, "t0");
    assert_3D_vector(t1, "t1");
    assert_3D_vector(t2, "t2");

    Eigen::Vector3<npe_Scalar_p> p_copy, t0_copy, t1_copy, t2_copy;
    copy_vector(p, p_copy);
    copy_vector(t0, t0_copy);
    copy_vector(t1, t1_copy);
    copy_vector(t2, t2_copy);

    Eigen::Matrix<npe_Scalar_p, 12, 12> hess;
    ipc::point_plane_distance_hessian(p_copy, t0_copy, t1_copy, t2_copy, hess);
    return npe::move(hess);

npe_end_code()
