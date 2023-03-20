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
.. doxygenfunction:: edge_edge_mollifier(const Eigen::Ref<const Eigen::Vector3d>& ea0, const Eigen::Ref<const Eigen::Vector3d>& ea1, const Eigen::Ref<const Eigen::Vector3d>& eb0, const Eigen::Ref<const Eigen::Vector3d>& eb1, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier(const double x, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier_gradient(const Eigen::Ref<const Eigen::Vector3d>& ea0, const Eigen::Ref<const Eigen::Vector3d>& ea1, const Eigen::Ref<const Eigen::Vector3d>& eb0, const Eigen::Ref<const Eigen::Vector3d>& eb1, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier_gradient(const double x, const double eps_x)
.. doxygenfunction:: edge_edge_mollifier_hessian(const Eigen::Ref<const Eigen::Vector3d>& ea0, const Eigen::Ref<const Eigen::Vector3d>& ea1, const Eigen::Ref<const Eigen::Vector3d>& eb0, const Eigen::Ref<const Eigen::Vector3d>& eb1, const double eps_x)
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

.. doxygenfunction:: point_plane_distance(const Eigen::Ref<const Eigen::Vector3d>& p, const Eigen::Ref<const Eigen::Vector3d>& origin, const Eigen::Ref<const Eigen::Vector3d>& normal)
.. doxygenfunction:: point_plane_distance(const Eigen::Ref<const Eigen::Vector3d>& p, const Eigen::Ref<const Eigen::Vector3d>& t0, const Eigen::Ref<const Eigen::Vector3d>& t1, const Eigen::Ref<const Eigen::Vector3d>& t2)
.. doxygenfunction:: point_plane_distance_gradient(const Eigen::Ref<const Eigen::Vector3d>& p, const Eigen::Ref<const Eigen::Vector3d>& origin, const Eigen::Ref<const Eigen::Vector3d>& normal)
.. doxygenfunction:: point_plane_distance_gradient(const Eigen::Ref<const Eigen::Vector3d>& p, const Eigen::Ref<const Eigen::Vector3d>& t0, const Eigen::Ref<const Eigen::Vector3d>& t1, const Eigen::Ref<const Eigen::Vector3d>& t2)
.. doxygenfunction:: point_plane_distance_hessian(const Eigen::Ref<const Eigen::Vector3d>& p, const Eigen::Ref<const Eigen::Vector3d>& origin, const Eigen::Ref<const Eigen::Vector3d>& normal)
.. doxygenfunction:: point_plane_distance_hessian(const Eigen::Ref<const Eigen::Vector3d>& p, const Eigen::Ref<const Eigen::Vector3d>& t0, const Eigen::Ref<const Eigen::Vector3d>& t1, const Eigen::Ref<const Eigen::Vector3d>& t2)

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