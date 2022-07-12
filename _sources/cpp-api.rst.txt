API
===

Main Functions
--------------
.. doxygenfunction:: construct_constraint_set(const CollisionMesh& mesh, const Eigen::MatrixXd& V, const double dhat, Constraints& constraint_set, const double dmin = 0, const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID)
.. doxygenfunction:: construct_constraint_set(const Candidates& candidates, const CollisionMesh& mesh, const Eigen::MatrixXd& V, const double dhat, Constraints& constraint_set, const double dmin = 0)

.. doxygenfunction:: compute_barrier_potential
.. doxygenfunction:: compute_barrier_potential_gradient
.. doxygenfunction:: compute_barrier_potential_hessian

.. doxygenfunction:: is_step_collision_free(const CollisionMesh& mesh,const Eigen::MatrixXd& V0,const Eigen::MatrixXd& V1,const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID,const double tolerance = 1e-6,const long max_iterations = 1e7)
.. doxygenfunction:: is_step_collision_free(const Candidates& candidates,const CollisionMesh& mesh,const Eigen::MatrixXd& V0,const Eigen::MatrixXd& V1,const double tolerance = 1e-6,const long max_iterations = 1e7)

.. doxygenfunction:: compute_collision_free_stepsize(const CollisionMesh& mesh,const Eigen::MatrixXd& V0,const Eigen::MatrixXd& V1,const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID,const double tolerance = 1e-6,const long max_iterations = 1e7)
.. doxygenfunction:: compute_collision_free_stepsize(const Candidates& candidates,const CollisionMesh& mesh,const Eigen::MatrixXd& V0,const Eigen::MatrixXd& V1,const double tolerance = 1e-6,const long max_iterations = 1e7)

Collision Mesh
--------------

.. doxygenclass:: ipc::CollisionMesh

Collision Constraints
---------------------

.. doxygenstruct:: ipc::CollisionConstraint
.. doxygenstruct:: ipc::VertexVertexConstraint
.. doxygenstruct:: ipc::EdgeVertexConstraint
.. doxygenstruct:: ipc::EdgeEdgeConstraint
.. doxygenstruct:: ipc::FaceVertexConstraint
.. doxygenstruct:: ipc::PlaneVertexConstraint
.. doxygenstruct:: ipc::Constraints

Barrier
-------

.. doxygenfunction:: barrier

.. doxygenfunction:: barrier_gradient

.. doxygenfunction:: barrier_hessian

Distance
--------

Distance Type
^^^^^^^^^^^^^

.. doxygenenum:: PointEdgeDistanceType
.. doxygenenum:: EdgeEdgeDistanceType
.. doxygenenum:: PointTriangleDistanceType

.. doxygenfunction:: point_edge_distance_type
.. doxygenfunction:: edge_edge_distance_type
.. doxygenfunction:: point_triangle_distance_type

Edge-Edge Mollifier
^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^

.. doxygenfunction:: edge_edge_distance(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1)
.. doxygenfunction:: edge_edge_distance(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const EdgeEdgeDistanceType dtype)
.. doxygenfunction:: edge_edge_distance_gradient(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: edge_edge_distance_gradient(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const EdgeEdgeDistanceType dtype, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: edge_edge_distance_hessian(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: edge_edge_distance_hessian(const Eigen::MatrixBase<DerivedEA0> &ea0, const Eigen::MatrixBase<DerivedEA1> &ea1, const Eigen::MatrixBase<DerivedEB0> &eb0, const Eigen::MatrixBase<DerivedEB1> &eb1, const EdgeEdgeDistanceType dtype, Eigen::PlainObjectBase<DerivedHess> &hess)

Line-Line
^^^^^^^^^

.. doxygenfunction:: line_line_distance
.. doxygenfunction:: ipc::line_line_distance_gradient
.. doxygenfunction:: ipc::line_line_distance_hessian

Point-Edge
^^^^^^^^^^

.. doxygenfunction:: point_edge_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1)
.. doxygenfunction:: point_edge_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, const PointEdgeDistanceType dtype)
.. doxygenfunction:: point_edge_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_edge_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, const PointEdgeDistanceType dtype, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_edge_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: point_edge_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedE0> &e0, const Eigen::MatrixBase<DerivedE1> &e1, const PointEdgeDistanceType dtype, Eigen::PlainObjectBase<DerivedHess> &hess)

Point-Line
^^^^^^^^^^

.. doxygenfunction:: point_line_distance
.. doxygenfunction:: point_line_distance_gradient
.. doxygenfunction:: point_line_distance_hessian

Point-Plane
^^^^^^^^^^^

.. doxygenfunction:: point_plane_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedOrigin> &origin, const Eigen::MatrixBase<DerivedNormal> &normal)
.. doxygenfunction:: point_plane_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2)
.. doxygenfunction:: point_plane_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedOrigin> &origin, const Eigen::MatrixBase<DerivedNormal> &normal, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_plane_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_plane_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedOrigin> &origin, const Eigen::MatrixBase<DerivedNormal> &normal, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: point_plane_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, Eigen::PlainObjectBase<DerivedHess> &hess)

Point-Point
^^^^^^^^^^^

.. doxygenfunction:: point_point_distance
.. doxygenfunction:: point_point_distance_gradient
.. doxygenfunction:: point_point_distance_hessian

Point-Triangle
^^^^^^^^^^^^^^

.. doxygenfunction:: point_triangle_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2)
.. doxygenfunction:: point_triangle_distance(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, const PointTriangleDistanceType dtype)
.. doxygenfunction:: point_triangle_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_triangle_distance_gradient(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, const PointTriangleDistanceType dtype, Eigen::PlainObjectBase<DerivedGrad> &grad)
.. doxygenfunction:: point_triangle_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, Eigen::PlainObjectBase<DerivedHess> &hess)
.. doxygenfunction:: point_triangle_distance_hessian(const Eigen::MatrixBase<DerivedP> &p, const Eigen::MatrixBase<DerivedT0> &t0, const Eigen::MatrixBase<DerivedT1> &t1, const Eigen::MatrixBase<DerivedT2> &t2, const PointTriangleDistanceType dtype, Eigen::PlainObjectBase<DerivedHess> &hess)

CCD
---

Broad-Phase
^^^^^^^^^^^

.. doxygenenum:: ipc::BroadPhaseMethod

Candidates
""""""""""

.. doxygenstruct:: ipc::ContinuousCollisionCandidate

.. doxygenstruct:: ipc::VertexVertexCandidate
.. doxygenstruct:: ipc::EdgeVertexCandidate
.. doxygenstruct:: ipc::EdgeEdgeCandidate
.. doxygenstruct:: ipc::EdgeFaceCandidate
.. doxygenstruct:: ipc::FaceVertexCandidate

.. doxygenstruct:: ipc::Candidates

Narrow-Phase
^^^^^^^^^^^^

.. doxygenvariable:: ipc::DEFAULT_CCD_TOLERANCE
.. doxygenvariable:: ipc::DEFAULT_CCD_MAX_ITERATIONS
.. doxygenvariable:: ipc::DEFAULT_CCD_CONSERVATIVE_RESCALING

Utils
-----

.. doxygenfunction:: compute_minimum_distance

.. doxygenfunction:: has_intersections
