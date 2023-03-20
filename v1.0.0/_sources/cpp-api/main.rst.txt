Main Functions
==============

.. doxygenfunction:: construct_friction_constraint_set(const CollisionMesh& mesh, const Eigen::MatrixXd& V, const Constraints& contact_constraint_set, double dhat, double barrier_stiffness, double mu, FrictionConstraints& friction_constraint_set)
.. doxygenfunction:: construct_friction_constraint_set(const CollisionMesh& mesh, const Eigen::MatrixXd& V, const Constraints& contact_constraint_set, double dhat, double barrier_stiffness, const Eigen::VectorXd& mus, FrictionConstraints& friction_constraint_set)
.. doxygenfunction:: construct_friction_constraint_set( const CollisionMesh& mesh, const Eigen::MatrixXd& V, const Constraints& contact_constraint_set, double dhat, double barrier_stiffness, const Eigen::VectorXd& mus, const std::function<double(double, double)>& blend_mu, FrictionConstraints& friction_constraint_set)

.. doxygenfunction:: compute_barrier_potential
.. doxygenfunction:: compute_barrier_potential_gradient
.. doxygenfunction:: compute_barrier_potential_hessian

.. doxygenfunction:: compute_friction_potential
.. doxygenfunction:: compute_friction_potential_gradient
.. doxygenfunction:: compute_friction_potential_hessian

.. doxygenfunction:: is_step_collision_free(const CollisionMesh& mesh,const Eigen::MatrixXd& V0,const Eigen::MatrixXd& V1,const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID,const double tolerance = 1e-6,const long max_iterations = 1e7)
.. doxygenfunction:: is_step_collision_free(const Candidates& candidates,const CollisionMesh& mesh,const Eigen::MatrixXd& V0,const Eigen::MatrixXd& V1,const double tolerance = 1e-6,const long max_iterations = 1e7)

.. doxygenfunction:: compute_collision_free_stepsize(const CollisionMesh& mesh,const Eigen::MatrixXd& V0,const Eigen::MatrixXd& V1,const BroadPhaseMethod method = BroadPhaseMethod::HASH_GRID,const double tolerance = 1e-6,const long max_iterations = 1e7)
.. doxygenfunction:: compute_collision_free_stepsize(const Candidates& candidates,const CollisionMesh& mesh,const Eigen::MatrixXd& V0,const Eigen::MatrixXd& V1,const double tolerance = 1e-6,const long max_iterations = 1e7)
