Distance
========

Distance Type
-------------

.. doxygenenum:: PointEdgeDistanceType
.. doxygenenum:: EdgeEdgeDistanceType
.. doxygenenum:: PointTriangleDistanceType

.. doxygenfunction:: point_edge_distance_type
.. doxygenfunction:: edge_edge_distance_type
.. doxygenfunction:: point_triangle_distance_type

Edge-Edge Mollifier
-------------------

.. doxygenfunction:: edge_edge_mollifier_threshold
.. doxygenfunction:: edge_edge_cross_squarednorm
.. doxygenfunction:: edge_edge_cross_squarednorm_gradient(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: edge_edge_cross_squarednorm_hessian(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: edge_edge_mollifier(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier(const T &x, double eps_x)
.. doxygenfunction:: edge_edge_mollifier_gradient(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const double eps_x, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: edge_edge_mollifier_gradient(const T &x, double eps_x)
.. doxygenfunction:: edge_edge_mollifier_hessian(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const double eps_x, Eigen::PlainObjectBase<DerivedHess> &hess))
.. doxygenfunction:: edge_edge_mollifier_hessian(const T &x, double eps_x)

Edge-Edge
---------

.. doxygenfunction:: edge_edge_distance(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1)
.. doxygenfunction:: edge_edge_distance(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const EdgeEdgeDistanceType dtype)
.. doxygenfunction:: edge_edge_distance_gradient(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: edge_edge_distance_gradient(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const EdgeEdgeDistanceType dtype, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: edge_edge_distance_hessian(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: edge_edge_distance_hessian(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const EdgeEdgeDistanceType dtype, Eigen::PlainObjectBase<DerivedHess> &hess)

Line-Line
---------

.. doxygenfunction:: line_line_distance
.. doxygenfunction:: ipc::line_line_distance_gradient
.. doxygenfunction:: ipc::line_line_distance_hessian

Point-Edge
----------

.. doxygenfunction:: point_edge_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1)
.. doxygenfunction:: point_edge_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, const PointEdgeDistanceType dtype)
.. doxygenfunction:: point_edge_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_edge_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, const PointEdgeDistanceType dtype, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_edge_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: point_edge_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, const PointEdgeDistanceType dtype, Eigen::PlainObjectBase<DerivedHess> &hess)

Point-Line
----------

.. doxygenfunction:: point_line_distance
.. doxygenfunction:: point_line_distance_gradient
.. doxygenfunction:: point_line_distance_hessian

Point-Plane
-----------

.. doxygenfunction:: point_plane_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedOrigin> &origin, const Eigen::MatrixBase<DerivedNormal> &normal)
.. doxygenfunction:: point_plane_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2)
.. doxygenfunction:: point_plane_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedOrigin> &origin, const Eigen::MatrixBase<DerivedNormal> &normal, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_plane_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_plane_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedOrigin> &origin, const Eigen::MatrixBase<DerivedNormal> &normal, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: point_plane_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, Eigen::PlainObjectBase<DerivedHess> &hess)

Point-Point
-----------

.. doxygenfunction:: point_point_distance
.. doxygenfunction:: point_point_distance_gradient
.. doxygenfunction:: point_point_distance_hessian

Point-Triangle
--------------

.. doxygenfunction:: point_triangle_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2)
.. doxygenfunction:: point_triangle_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, const PointTriangleDistanceType dtype)
.. doxygenfunction:: point_triangle_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_triangle_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, const PointTriangleDistanceType dtype, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_triangle_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: point_triangle_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, const PointTriangleDistanceType dtype, Eigen::PlainObjectBase<DerivedHess> &hess)