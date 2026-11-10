I/O
===

glTF
----

Reading and writing glTF scenes requires building with the
``IPC_TOOLKIT_WITH_GLTF`` CMake option.

.. doxygenfunction:: ipc::rigid::read_gltf
.. doxygenfunction:: ipc::rigid::write_gltf

OBJ
---

.. doxygenfunction:: ipc::write_candidates_obj(std::ostream &out, Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Eigen::MatrixXi> E, Eigen::ConstRef<Eigen::MatrixXi> F, const CandidateVector<Candidate> &candidates, const int v_offset)
.. doxygenfunction:: ipc::write_candidates_obj(const std::string &filename, Eigen::ConstRef<Eigen::MatrixXd> V, Eigen::ConstRef<Eigen::MatrixXi> E, Eigen::ConstRef<Eigen::MatrixXi> F, const CandidateVector<Candidate> &candidates)
