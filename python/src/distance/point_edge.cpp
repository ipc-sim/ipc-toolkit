// clang-format off

#include <npe.h>

#include <ipc/distance/point_edge.hpp>

#include <common.hpp>

const char* ds_point_edge_distance = R"ipc_Qu8mg5v7(
Compute the distance between a point and edge in 2D or 3D.

Parameters
----------
p : The point
e0 : The first vertex of the edge
e1 : The second vertex of the edge
dtype : (Optional) The point edge distance type to compute

Returns
-------
The distance between the point and edge

See also
--------
point_edge_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_edge_distance)
npe_doc(ds_point_edge_distance)
npe_arg(p, dense_double)
npe_arg(e0, npe_matches(p))
npe_arg(e1, npe_matches(p))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::PointEdgeDistanceType::P_E0) == 0, "PointEdgeDistanceType enum changed!");
    static_assert(int(ipc::PointEdgeDistanceType::P_E1) == 1, "PointEdgeDistanceType enum changed!");
    static_assert(int(ipc::PointEdgeDistanceType::P_E) == 2, "PointEdgeDistanceType enum changed!");

    assert_2D_or_3D_vector(p, "p");
    assert_2D_or_3D_vector(e0, "e0");
    assert_2D_or_3D_vector(e1, "e1");

    Eigen::VectorX3<npe_Scalar_p> p_copy, e0_copy, e1_copy;
    copy_vector(p, p_copy);
    copy_vector(e0, e0_copy);
    copy_vector(e1, e1_copy);

    if(dtype < 0) {
        return ipc::point_edge_distance(p_copy, e0_copy, e1_copy);
    } else if(dtype < 2) {
        return ipc::point_edge_distance(p_copy, e0_copy, e1_copy, ipc::PointEdgeDistanceType(dtype));
    } else {
        std::stringstream ss;
        ss << "invalid point-edge dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_edge_distance_gradient = R"ipc_Qu8mg5v7(
Compute the gradient of the distance between a point and edge.

Parameters
----------
p : The point
e0 : The first vertex of the edge
e1 : The second vertex of the edge
dtype : (Optional) The point edge distance type to compute

Returns
-------
The gradient of the distance wrt p, e0, and e1.

See also
--------
point_edge_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_edge_distance_gradient)
npe_doc(ds_point_edge_distance_gradient)
npe_arg(p, dense_double)
npe_arg(e0, npe_matches(p))
npe_arg(e1, npe_matches(p))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::PointEdgeDistanceType::P_E0) == 0, "PointEdgeDistanceType enum changed!");
    static_assert(int(ipc::PointEdgeDistanceType::P_E1) == 1, "PointEdgeDistanceType enum changed!");
    static_assert(int(ipc::PointEdgeDistanceType::P_E) == 2, "PointEdgeDistanceType enum changed!");

    assert_2D_or_3D_vector(p, "p");
    assert_2D_or_3D_vector(e0, "e0");
    assert_2D_or_3D_vector(e1, "e1");

    Eigen::VectorX3<npe_Scalar_p> p_copy, e0_copy, e1_copy;
    copy_vector(p, p_copy);
    copy_vector(e0, e0_copy);
    copy_vector(e1, e1_copy);

    Eigen::VectorX12<npe_Scalar_p> grad;
    if(dtype < 0) {
        ipc::point_edge_distance_gradient(p_copy, e0_copy, e1_copy, grad);
    } else if(dtype < 2) {
        ipc::point_edge_distance_gradient(p_copy, e0_copy, e1_copy, ipc::PointEdgeDistanceType(dtype), grad);
    } else {
        std::stringstream ss;
        ss << "invalid point-edge dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }
    return npe::move(grad);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_edge_distance_hessian = R"ipc_Qu8mg5v7(
Compute the hessian of the distance between a point and edge.

Parameters
----------
p : The point
e0 : The first vertex of the edge
e1 : The second vertex of the edge
dtype : (Optional) The point edge distance type to compute

Returns
-------
The hessian of the distance wrt p, e0, and e1.

See also
--------
point_edge_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_edge_distance_hessian)
npe_doc(ds_point_edge_distance_hessian)
npe_arg(p, dense_double)
npe_arg(e0, npe_matches(p))
npe_arg(e1, npe_matches(p))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::PointEdgeDistanceType::P_E0) == 0, "PointEdgeDistanceType enum changed!");
    static_assert(int(ipc::PointEdgeDistanceType::P_E1) == 1, "PointEdgeDistanceType enum changed!");
    static_assert(int(ipc::PointEdgeDistanceType::P_E) == 2, "PointEdgeDistanceType enum changed!");

    assert_2D_or_3D_vector(p, "p");
    assert_2D_or_3D_vector(e0, "e0");
    assert_2D_or_3D_vector(e1, "e1");

    Eigen::VectorX3<npe_Scalar_p> p_copy, e0_copy, e1_copy;
    copy_vector(p, p_copy);
    copy_vector(e0, e0_copy);
    copy_vector(e1, e1_copy);

    Eigen::MatrixXX12<npe_Scalar_p> hess;
    if(dtype < 0) {
        ipc::point_edge_distance_hessian(p_copy, e0_copy, e1_copy, hess);
    } else if(dtype < 2) {
        ipc::point_edge_distance_hessian(p_copy, e0_copy, e1_copy, ipc::PointEdgeDistanceType(dtype), hess);
    } else {
        std::stringstream ss;
        ss << "invalid point-edge dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }
    return npe::move(hess);

npe_end_code()
