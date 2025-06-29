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
.. doxygenfunction:: ipc::edge_edge_cross_squarednorm_gradient
.. doxygenfunction:: ipc::edge_edge_cross_squarednorm_hessian
.. doxygenfunction:: edge_edge_mollifier(Eigen::ConstRef<Eigen::Vector3d> ea0, Eigen::ConstRef<Eigen::Vector3d> ea1, Eigen::ConstRef<Eigen::Vector3d> eb0, Eigen::ConstRef<Eigen::Vector3d> eb1, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier(const double x, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier_gradient(Eigen::ConstRef<Eigen::Vector3d> ea0, Eigen::ConstRef<Eigen::Vector3d> ea1, Eigen::ConstRef<Eigen::Vector3d> eb0, Eigen::ConstRef<Eigen::Vector3d> eb1, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier_gradient(const double x, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier_hessian(Eigen::ConstRef<Eigen::Vector3d> ea0, Eigen::ConstRef<Eigen::Vector3d> ea1, Eigen::ConstRef<Eigen::Vector3d> eb0, Eigen::ConstRef<Eigen::Vector3d> eb1, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier_hessian(const double x, const double eps_x)

Edge-Edge
---------

.. doxygenfunction:: edge_edge_distance
.. doxygenfunction:: edge_edge_distance_gradient
.. doxygenfunction:: edge_edge_distance_hessian

Line-Line
---------

.. doxygenfunction:: line_line_distance
.. doxygenfunction:: ipc::line_line_distance_gradient
.. doxygenfunction:: ipc::line_line_distance_hessian

Point-Edge
----------

.. doxygenfunction:: point_edge_distance
.. doxygenfunction:: point_edge_distance_gradient
.. doxygenfunction:: point_edge_distance_hessian

Point-Line
----------

.. doxygenfunction:: point_line_distance
.. doxygenfunction:: point_line_distance_gradient
.. doxygenfunction:: point_line_distance_hessian

Point-Plane
-----------

.. doxygenfunction:: point_plane_distance(Eigen::ConstRef<Eigen::Vector3d> p, Eigen::ConstRef<Eigen::Vector3d> origin, Eigen::ConstRef<Eigen::Vector3d> normal)
.. doxygenfunction:: point_plane_distance(Eigen::ConstRef<Eigen::Vector3d> p, Eigen::ConstRef<Eigen::Vector3d> t0, Eigen::ConstRef<Eigen::Vector3d> t1, Eigen::ConstRef<Eigen::Vector3d> t2)
.. doxygenfunction:: point_plane_distance_gradient(Eigen::ConstRef<Eigen::Vector3d> p, Eigen::ConstRef<Eigen::Vector3d> origin, Eigen::ConstRef<Eigen::Vector3d> normal)
.. doxygenfunction:: point_plane_distance_gradient(Eigen::ConstRef<Eigen::Vector3d> p, Eigen::ConstRef<Eigen::Vector3d> t0, Eigen::ConstRef<Eigen::Vector3d> t1, Eigen::ConstRef<Eigen::Vector3d> t2)
.. doxygenfunction:: point_plane_distance_hessian(Eigen::ConstRef<Eigen::Vector3d> p, Eigen::ConstRef<Eigen::Vector3d> origin, Eigen::ConstRef<Eigen::Vector3d> normal)
.. doxygenfunction:: point_plane_distance_hessian(Eigen::ConstRef<Eigen::Vector3d> p, Eigen::ConstRef<Eigen::Vector3d> t0, Eigen::ConstRef<Eigen::Vector3d> t1, Eigen::ConstRef<Eigen::Vector3d> t2)

Point-Point
-----------

.. doxygenfunction:: point_point_distance
.. doxygenfunction:: point_point_distance_gradient
.. doxygenfunction:: point_point_distance_hessian

Point-Triangle
--------------

.. doxygenfunction:: point_triangle_distance
.. doxygenfunction:: point_triangle_distance_gradient
.. doxygenfunction:: point_triangle_distance_hessian