// clang-format off

#include <npe.h>

#include <ipc/distance/edge_edge.hpp>

#include <common.hpp>

const char* ds_edge_edge_distance = R"ipc_Qu8mg5v7(
Compute the distance between a two lines segments in 3D.

Parameters
----------
ea0 : The first vertex of the first edge
ea1 : The second vertex of the first edge
eb0 : The first vertex of the second edge
eb1 : The second vertex of the second edge
dtype : (Optional) The edge-edge distance type to compute

Returns
-------
The distance between the two edges.

See also
--------
edge_edge_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_distance)
npe_doc(ds_edge_edge_distance)
npe_arg(ea0, dense_double)
npe_arg(ea1, npe_matches(ea0))
npe_arg(eb0, npe_matches(ea0))
npe_arg(eb1, npe_matches(ea0))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB0) == 0, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB1) == 1, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB0) == 2, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB1) == 3, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB0) == 4, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB1) == 5, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB) == 6, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB) == 7, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB) == 8, "EdgeEdgeDistanceType enum changed!");

    assert_3D_vector(ea0, "ea0");
    assert_3D_vector(ea1, "ea1");
    assert_3D_vector(eb0, "eb0");
    assert_3D_vector(eb1, "eb1");

    Eigen::Vector3<npe_Scalar_ea0> ea0_copy, ea1_copy, eb0_copy, eb1_copy;
    copy_vector(ea0, ea0_copy);
    copy_vector(ea1, ea1_copy);
    copy_vector(eb0, eb0_copy);
    copy_vector(eb1, eb1_copy);

    if(dtype < 0) {
        return ipc::edge_edge_distance(ea0_copy, ea1_copy, eb0_copy, eb1_copy);
    } else if(dtype < 2) {
        return ipc::edge_edge_distance(
            ea0_copy, ea1_copy, eb0_copy, eb1_copy,
            ipc::EdgeEdgeDistanceType(dtype));
    } else {
        std::stringstream ss;
        ss << "invalid edge-edge dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_distance_gradient = R"ipc_Qu8mg5v7(
Compute the gradient of the distance between a two lines segments.

Parameters
----------
ea0 : The first vertex of the first edge
ea1 : The second vertex of the first edge
eb0 : The first vertex of the second edge
eb1 : The second vertex of the second edge
dtype : (Optional) The point edge distance type to compute

Returns
-------
The gradient of the distance wrt ea0, ea1, eb0, and eb1.

See also
--------
edge_edge_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_distance_gradient)
npe_doc(ds_edge_edge_distance_gradient)
npe_arg(ea0, dense_double)
npe_arg(ea1, npe_matches(ea0))
npe_arg(eb0, npe_matches(ea0))
npe_arg(eb1, npe_matches(ea0))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB0) == 0, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB1) == 1, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB0) == 2, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB1) == 3, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB0) == 4, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB1) == 5, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB) == 6, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB) == 7, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB) == 8, "EdgeEdgeDistanceType enum changed!");

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
    if(dtype < 0) {
        ipc::edge_edge_distance_gradient(
            ea0_copy, ea1_copy, eb0_copy, eb1_copy, grad);
    } else if(dtype < 2) {
        ipc::edge_edge_distance_gradient(
            ea0_copy, ea1_copy, eb0_copy, eb1_copy,
            ipc::EdgeEdgeDistanceType(dtype), grad);
    } else {
        std::stringstream ss;
        ss << "invalid edge-edge dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }
    return npe::move(grad);

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_distance_hessian = R"ipc_Qu8mg5v7(
Compute the hessian of the distance between a two lines segments.

Parameters
----------
ea0 : The first vertex of the first edge
ea1 : The second vertex of the first edge
eb0 : The first vertex of the second edge
eb1 : The second vertex of the second edge
dtype : (Optional) The point edge distance type to compute

Returns
-------
The hessian of the distance wrt ea0, ea1, eb0, and eb1.

See also
--------
edge_edge_distance_type

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_distance_hessian)
npe_doc(ds_edge_edge_distance_hessian)
npe_arg(ea0, dense_double)
npe_arg(ea1, npe_matches(ea0))
npe_arg(eb0, npe_matches(ea0))
npe_arg(eb1, npe_matches(ea0))
npe_default_arg(dtype, int, -1)
npe_begin_code()

    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB0) == 0, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB1) == 1, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB0) == 2, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB1) == 3, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB0) == 4, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB1) == 5, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA0_EB) == 6, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA1_EB) == 7, "EdgeEdgeDistanceType enum changed!");
    static_assert(int(ipc::EdgeEdgeDistanceType::EA_EB) == 8, "EdgeEdgeDistanceType enum changed!");

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
    if(dtype < 0) {
        ipc::edge_edge_distance_hessian(
            ea0_copy, ea1_copy, eb0_copy, eb1_copy, hess);
    } else if(dtype < 2) {
        ipc::edge_edge_distance_hessian(
            ea0_copy, ea1_copy, eb0_copy, eb1_copy,
            ipc::EdgeEdgeDistanceType(dtype), hess);
    } else {
        std::stringstream ss;
        ss << "invalid edge-edge dtype (" << dtype << ")";
        throw pybind11::value_error(ss.str());
    }
    return npe::move(hess);

npe_end_code()
