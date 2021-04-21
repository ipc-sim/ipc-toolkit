// clang-format off

#include <npe.h>

#include <ipc/distance/point_edge.hpp>

#include <common.hpp>

const char* ds_point_edge_distance_type = R"ipc_Qu8mg5v7(
Determine the closest pair between a point and edge.

Parameters
----------
p  : The point
e0 : The first vertex of the edge
e1 : The second vertex of the edge

Returns
-------
The distance type of the point-edge pair.

0: P_E0
1: P_E1
2: P_E

See also
--------

Notes
-----

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_edge_distance_type)
npe_doc(ds_point_edge_distance_type)
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

    return int(ipc::point_edge_distance_type(p_copy, e0_copy, e1_copy));

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_point_triangle_distance_type = R"ipc_Qu8mg5v7(
Determine the closest pair between a point and triangle.

Parameters
----------
p  : The point
t0 : The first vertex of the triangle
t1 : The second vertex of the triangle
t2 : The third vertex of the triangle

Returns
-------
The distance type of the point-triangle pair.

0: P_T0
1: P_T1
2: P_T2
3: P_E0
4: P_E1
5: P_E2
6: P_T

See also
--------

Notes
-----

Examples
--------

)ipc_Qu8mg5v7";

npe_function(point_triangle_distance_type)
npe_doc(ds_point_triangle_distance_type)
npe_arg(p, dense_double)
npe_arg(t0, npe_matches(p))
npe_arg(t1, npe_matches(p))
npe_arg(t2, npe_matches(p))
npe_begin_code()

    assert_3D_vector(p, "p");
    assert_3D_vector(t0, "t0");
    assert_3D_vector(t1, "t1");
    assert_3D_vector(t2, "t2");

    Eigen::VectorX3<npe_Scalar_p> p_copy, t0_copy, t1_copy, t2_copy;
    copy_vector(p, p_copy);
    copy_vector(t0, t0_copy);
    copy_vector(t1, t1_copy);
    copy_vector(t2, t2_copy);

    return int(ipc::point_triangle_distance_type(p_copy, t0_copy, t1_copy, t2_copy));

npe_end_code()

///////////////////////////////////////////////////////////////////////////////

const char* ds_edge_edge_distance_type = R"ipc_Qu8mg5v7(
Determine the closest pair between two edges.

Parameters
----------
ea0 : The first vertex of the first edge
ea1 : The second vertex of the first edge
eb0 : The first vertex of the second edge
eb1 : The second vertex of the second edge

Returns
-------
The distance type of the edge-edge pair.

0: EA0_EB0
1: EA0_EB1
2: EA1_EB0
3: EA1_EB1
4: EA_EB0
5: EA_EB1
6: EA0_EB
7: EA1_EB
8: EA_EB

See also
--------

Notes
-----
The distance is actually squared distance.

Examples
--------

)ipc_Qu8mg5v7";

npe_function(edge_edge_distance_type)
npe_doc(ds_edge_edge_distance_type)
npe_arg(ea0, dense_double)
npe_arg(ea1, npe_matches(ea0))
npe_arg(eb0, npe_matches(ea0))
npe_arg(eb1, npe_matches(ea0))
npe_begin_code()

    assert_3D_vector(ea0, "ea0");
    assert_3D_vector(ea1, "ea1");
    assert_3D_vector(eb0, "eb0");
    assert_3D_vector(eb1, "eb1");

    Eigen::VectorX3<npe_Scalar_ea0> ea0_copy, ea1_copy, eb0_copy, eb1_copy;
    copy_vector(ea0, ea0_copy);
    copy_vector(ea1, ea1_copy);
    copy_vector(eb0, eb0_copy);
    copy_vector(eb1, eb1_copy);

    return int(ipc::point_triangle_distance_type(ea0_copy, ea1_copy, eb0_copy, eb1_copy));

npe_end_code()
