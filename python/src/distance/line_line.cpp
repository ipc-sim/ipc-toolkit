// clang-format off

#include <npe.h>

#include <ipc/distance/line_line.hpp>

#include <common.hpp>

const char* ds_line_line_distance = R"ipc_Qu8mg5v7(
Compute the distance between a two infinite lines in 3D.

Parameters
----------
ea0 : The first vertex of the edge defining the first line
ea1 : The second vertex of the edge defining the first line
eb0 : The first vertex of the edge defining the second line
eb1 : The second vertex of the edge defining the second line

Returns
-------
The distance between the two lines.

See also
--------

Notes
-----
The distance is actually squared distance.

Warning
-------
If the lines are parallel this function returns a distance of zero.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(line_line_distance)
npe_doc(ds_line_line_distance)
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

    return ipc::line_line_distance(ea0_copy, ea1_copy, eb0_copy, eb1_copy);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_line_line_distance_gradient = R"ipc_Qu8mg5v7(
Compute the gradient of the distance between a two lines in 3D.

Parameters
----------
ea0 : The first vertex of the edge defining the first line
ea1 : The second vertex of the edge defining the first line
eb0 : The first vertex of the edge defining the second line
eb1 : The second vertex of the edge defining the second line

Returns
-------
The gradient of the distance wrt ea0, ea1, eb0, and eb1.

See also
--------

Notes
-----
The distance is actually squared distance.

Warning
-------
If the lines are parallel this function returns a distance of zero.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(line_line_distance_gradient)
npe_doc(ds_line_line_distance_gradient)
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
    ipc::line_line_distance_gradient(ea0_copy, ea1_copy, eb0_copy, eb1_copy, grad);
    return npe::move(grad);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_line_line_distance_hessian = R"ipc_Qu8mg5v7(
Compute the hessian of the distance between a two lines in 3D.

Parameters
----------
ea0 : The first vertex of the edge defining the first line
ea1 : The second vertex of the edge defining the first line
eb0 : The first vertex of the edge defining the second line
eb1 : The second vertex of the edge defining the second line

Returns
-------
The hessian of the distance wrt ea0, ea1, eb0, and eb1.

See also
--------

Notes
-----
The distance is actually squared distance.

Warning
-------
If the lines are parallel this function returns a distance of zero.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(line_line_distance_hessian)
npe_doc(ds_line_line_distance_hessian)
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
    ipc::line_line_distance_hessian(ea0_copy, ea1_copy, eb0_copy, eb1_copy, hess);
    return npe::move(hess);

npe_end_code()
