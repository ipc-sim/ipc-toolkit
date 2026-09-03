Distance
========

Distance Type
-------------

.. doxygenenum:: PointEdgeDistanceType
.. doxygenenum:: EdgeEdgeDistanceType
.. doxygenenum:: PointTriangleDistanceType

.. doxygenfunction:: ipc::point_edge_distance_type
.. doxygenfunction:: ipc::edge_edge_distance_type
.. doxygenfunction:: ipc::point_triangle_distance_type

Edge-Edge Mollifier
-------------------

.. doxygenfunction:: ipc::edge_edge_mollifier_threshold
.. doxygenfunction:: ipc::edge_edge_cross_squarednorm
.. doxygenfunction:: ipc::edge_edge_cross_squarednorm_gradient
.. doxygenfunction:: ipc::edge_edge_cross_squarednorm_hessian
.. doxygenfunction:: ipc::edge_edge_mollifier(const Eigen::MatrixBase<DerivedEA0>& ea0, const Eigen::MatrixBase<DerivedEA1>& ea1, const Eigen::MatrixBase<DerivedEB0>& eb0, const Eigen::MatrixBase<DerivedEB1>& eb1, const typename DerivedEA0::Scalar eps_x)
.. doxygenfunction:: ipc::edge_edge_mollifier(const T x, const T eps_x)
.. doxygenfunction:: ipc::edge_edge_mollifier_gradient(const Eigen::MatrixBase<DerivedEA0>& ea0, const Eigen::MatrixBase<DerivedEA1>& ea1, const Eigen::MatrixBase<DerivedEB0>& eb0, const Eigen::MatrixBase<DerivedEB1>& eb1, const typename DerivedEA0::Scalar eps_x)
.. doxygenfunction:: ipc::edge_edge_mollifier_gradient(const T x, const T eps_x)
.. doxygenfunction:: ipc::edge_edge_mollifier_hessian(const Eigen::MatrixBase<DerivedEA0>& ea0, const Eigen::MatrixBase<DerivedEA1>& ea1, const Eigen::MatrixBase<DerivedEB0>& eb0, const Eigen::MatrixBase<DerivedEB1>& eb1, const typename DerivedEA0::Scalar eps_x)
.. doxygenfunction:: ipc::edge_edge_mollifier_hessian(const T x, const T eps_x)

Edge-Edge
---------

.. doxygenfunction:: ipc::edge_edge_distance
.. doxygenfunction:: ipc::edge_edge_distance_gradient
.. doxygenfunction:: ipc::edge_edge_distance_hessian

Line-Line
---------

.. doxygenfunction:: ipc::line_line_distance
.. doxygenfunction:: ipc::line_line_distance_gradient
.. doxygenfunction:: ipc::line_line_distance_hessian

Point-Edge
----------

.. doxygenfunction:: ipc::point_edge_distance
.. doxygenfunction:: ipc::point_edge_distance_gradient
.. doxygenfunction:: ipc::point_edge_distance_hessian

Point-Line
----------

.. doxygenfunction:: ipc::point_line_distance
.. doxygenfunction:: ipc::point_line_distance_gradient
.. doxygenfunction:: ipc::point_line_distance_hessian

Point-Plane
-----------

.. doxygenfunction:: point_plane_distance(const Eigen::MatrixBase<DerivedP>& p, const Eigen::MatrixBase<DerivedOrigin>& origin, const Eigen::MatrixBase<DerivedNormal>& normal)
.. doxygenfunction:: point_plane_distance(const Eigen::MatrixBase<DerivedP>& p, const Eigen::MatrixBase<DerivedT0>& t0, const Eigen::MatrixBase<DerivedT1>& t1, const Eigen::MatrixBase<DerivedT2>& t2)
.. doxygenfunction:: point_plane_distance_gradient(const Eigen::MatrixBase<DerivedP>& p, const Eigen::MatrixBase<DerivedOrigin>& origin, const Eigen::MatrixBase<DerivedNormal>& normal)
.. doxygenfunction:: point_plane_distance_gradient(const Eigen::MatrixBase<DerivedP>& p, const Eigen::MatrixBase<DerivedT0>& t0, const Eigen::MatrixBase<DerivedT1>& t1, const Eigen::MatrixBase<DerivedT2>& t2)
.. doxygenfunction:: point_plane_distance_hessian(const Eigen::MatrixBase<DerivedP>& p, const Eigen::MatrixBase<DerivedOrigin>& origin, const Eigen::MatrixBase<DerivedNormal>& normal)
.. doxygenfunction:: point_plane_distance_hessian(const Eigen::MatrixBase<DerivedP>& p, const Eigen::MatrixBase<DerivedT0>& t0, const Eigen::MatrixBase<DerivedT1>& t1, const Eigen::MatrixBase<DerivedT2>& t2)

Point-Point
-----------

.. doxygenfunction:: ipc::point_point_distance
.. doxygenfunction:: ipc::point_point_distance_gradient
.. doxygenfunction:: ipc::point_point_distance_hessian

Point-Triangle
--------------

.. doxygenfunction:: ipc::point_triangle_distance
.. doxygenfunction:: ipc::point_triangle_distance_gradient
.. doxygenfunction:: ipc::point_triangle_distance_hessian