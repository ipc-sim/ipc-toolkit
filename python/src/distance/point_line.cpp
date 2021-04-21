// clang-format off

#include <npe.h>

#include <ipc/distance/point_line.hpp>

#include <common.hpp>

const char* ds_point_line_distance = R"ipc_Qu8mg5v7(
Compute the distance between a point and line in 2D or 3D.

Parameters
----------
p : The point
e0 : The first vertex of the edge defining the line
e1 : The second vertex of the edge defining the line

Returns
-------
The distance between the point and line.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_line_distance)
npe_doc(ds_point_line_distance)
npe_arg(p, dense_double)
npe_arg(e0, npe_matches(p))
npe_arg(e1, npe_matches(p))
npe_begin_code()

    assert_2D_or_3D_vector(p, "p");
    assert_2D_or_3D_vector(e0, "e0");
    assert_2D_or_3D_vector(e1, "e1");

    Eigen::VectorX3<npe_Scalar_p> p_copy, e0_copy, e1_copy;
    copy_vector(p, p_copy);
    copy_vector(e0, e0_copy);
    copy_vector(e1, e1_copy);

    return ipc::point_line_distance(p_copy, e0_copy, e1_copy);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_line_distance_gradient = R"ipc_Qu8mg5v7(
Compute the gradient of the distance between a point and line.

Parameters
----------
p : The point
e0 : The first vertex of the edge defining the line.
e1 : The second vertex of the edge defining the line.

Returns
-------
The gradient of the distance wrt p, e0, and e1.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_line_distance_gradient)
npe_doc(ds_point_line_distance_gradient)
npe_arg(p, dense_double)
npe_arg(e0, npe_matches(p))
npe_arg(e1, npe_matches(p))
npe_begin_code()

    assert_2D_or_3D_vector(p, "p");
    assert_2D_or_3D_vector(e0, "e0");
    assert_2D_or_3D_vector(e1, "e1");

    Eigen::VectorX3<npe_Scalar_p> p_copy, e0_copy, e1_copy;
    copy_vector(p, p_copy);
    copy_vector(e0, e0_copy);
    copy_vector(e1, e1_copy);

    Eigen::VectorX9<npe_Scalar_p> grad;
    ipc::point_line_distance_gradient(p_copy, e0_copy, e1_copy, grad);
    return npe::move(grad);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_line_distance_hessian = R"ipc_Qu8mg5v7(
Compute the hessian of the distance between a point and line.

Parameters
----------
p : The point
e0 : The first vertex of the edge defining the line
e1 : The second vertex of the edge defining the line

Returns
-------
The hessian of the distance wrt p, e0, and e1.

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_line_distance_hessian)
npe_doc(ds_point_line_distance_hessian)
npe_arg(p, dense_double)
npe_arg(e0, npe_matches(p))
npe_arg(e1, npe_matches(p))
npe_begin_code()

    assert_2D_or_3D_vector(p, "p");
    assert_2D_or_3D_vector(e0, "e0");
    assert_2D_or_3D_vector(e1, "e1");

    Eigen::VectorX3<npe_Scalar_p> p_copy, e0_copy, e1_copy;
    copy_vector(p, p_copy);
    copy_vector(e0, e0_copy);
    copy_vector(e1, e1_copy);

    Eigen::MatrixXX9<npe_Scalar_p> hess;
    ipc::point_line_distance_hessian(p_copy, e0_copy, e1_copy, hess);
    return npe::move(hess);

npe_end_code()
