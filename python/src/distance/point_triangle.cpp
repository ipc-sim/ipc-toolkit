// clang-format off

#include <npe.h>

#include <ipc/distance/point_triangle.hpp>

#include <common.hpp>

const char* ds_point_triangle_distance = R"ipc_Qu8mg5v7(
Compute the distance between a two lines segments in 3D.

Parameters
----------
p  : The point.
t0 : The first vertex of the triangle.
t1 : The second vertex of the triangle.
t2 : The third vertex of the triangle.
dtype : (Optional) The point-triangle distance type to compute

Returns
-------
The distance between the point and triangle.

See also
--------
point_triangle_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_triangle_distance)
npe_doc(ds_point_triangle_distance)
npe_arg(p, dense_double)
npe_arg(t0, npe_matches(p))
npe_arg(t1, npe_matches(p))
npe_arg(t2, npe_matches(p))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::PointTriangleDistanceType::P_T0) == 0, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T1) == 1, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T2) == 2, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E0) == 3, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E1) == 4, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E2) == 5, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T) == 6, "PointTriangleDistanceType enum changed!");

    assert_3D_vector(p, "p");
    assert_3D_vector(t0, "t0");
    assert_3D_vector(t1, "t1");
    assert_3D_vector(t2, "t2");

    Eigen::Vector3<npe_Scalar_p> p_copy, t0_copy, t1_copy, t2_copy;
    copy_vector(p, p_copy);
    copy_vector(t0, t0_copy);
    copy_vector(t1, t1_copy);
    copy_vector(t2, t2_copy);

    if(dtype < 0) {
        return ipc::point_triangle_distance(p_copy, t0_copy, t1_copy, t2_copy);
    } else if(dtype < 2) {
        return ipc::point_triangle_distance(
            p_copy, t0_copy, t1_copy, t2_copy,
            ipc::PointTriangleDistanceType(dtype));
    } else {
        std::stringstream ss;
        ss << "invalid point-triangle dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_triangle_distance_gradient = R"ipc_Qu8mg5v7(
Compute the gradient of the distance between a two lines segments.

Parameters
----------
p  : The point.
t0 : The first vertex of the triangle.
t1 : The second vertex of the triangle.
t2 : The third vertex of the triangle.
dtype : (Optional) The point-triangle distance type to compute

Returns
-------
The gradient of the distance wrt ea0, ea1, eb0, and eb1.

See also
--------
point_triangle_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_triangle_distance_gradient)
npe_doc(ds_point_triangle_distance_gradient)
npe_arg(p, dense_double)
npe_arg(t0, npe_matches(p))
npe_arg(t1, npe_matches(p))
npe_arg(t2, npe_matches(p))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::PointTriangleDistanceType::P_T0) == 0, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T1) == 1, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T2) == 2, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E0) == 3, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E1) == 4, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E2) == 5, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T) == 6, "PointTriangleDistanceType enum changed!");

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
    if(dtype < 0) {
        ipc::point_triangle_distance_gradient(
            p_copy, t0_copy, t1_copy, t2_copy, grad);
    } else if(dtype < 2) {
        ipc::point_triangle_distance_gradient(
            p_copy, t0_copy, t1_copy, t2_copy,
            ipc::PointTriangleDistanceType(dtype), grad);
    } else {
        std::stringstream ss;
        ss << "invalid point-triangle dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }
    return npe::move(grad);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_triangle_distance_hessian = R"ipc_Qu8mg5v7(
Compute the hessian of the distance between a two lines segments.

Parameters
----------
p  : The point.
t0 : The first vertex of the triangle.
t1 : The second vertex of the triangle.
t2 : The third vertex of the triangle.
dtype : (Optional) The point-triangle distance type to compute

Returns
-------
The hessian of the distance wrt ea0, ea1, eb0, and eb1.

See also
--------
point_triangle_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_triangle_distance_hessian)
npe_doc(ds_point_triangle_distance_hessian)
npe_arg(p, dense_double)
npe_arg(t0, npe_matches(p))
npe_arg(t1, npe_matches(p))
npe_arg(t2, npe_matches(p))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::PointTriangleDistanceType::P_T0) == 0, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T1) == 1, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T2) == 2, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E0) == 3, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E1) == 4, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_E2) == 5, "PointTriangleDistanceType enum changed!");
    static_assert(int(ipc::PointTriangleDistanceType::P_T) == 6, "PointTriangleDistanceType enum changed!");

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
    if(dtype < 0) {
        ipc::point_triangle_distance_hessian(
            p_copy, t0_copy, t1_copy, t2_copy, hess);
    } else if(dtype < 2) {
        ipc::point_triangle_distance_hessian(
            p_copy, t0_copy, t1_copy, t2_copy,
            ipc::PointTriangleDistanceType(dtype), hess);
    } else {
        std::stringstream ss;
        ss << "invalid point-triangle dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }
    return npe::move(hess);

npe_end_code()
