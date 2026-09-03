Utils
=====

Logger
------

.. doxygenfunction:: ipc::logger
.. doxygenfunction:: ipc::set_logger

Positive Semi-Definite Projection
---------------------------------

.. doxygenfunction:: ipc::project_to_psd
.. doxygenfunction:: ipc::project_to_pd

.. doxygenenum:: ipc::PSDProjectionMethod

Hessian Assembly
----------------

Pluggable backends for assembling per-collision Hessians into a global
matrix (see :cpp:func:`ipc::Potential::assemble_hessian`).

.. doxygenclass:: ipc::HessianAssembler
.. doxygenclass:: ipc::TripletHessianAssembler

The following backend is available when the toolkit is compiled with
``IPC_TOOLKIT_WITH_MESHFEM_SPARSE`` (the default), in which case it is also
what :cpp:func:`ipc::Potential::hessian` uses internally. It assembles into
`MeshFEMSparse <https://github.com/MeshFEM/MeshFEMSparse>`_'s block-CSC data
structures (no triplets, no ``setFromTriplets``) and reuses the sparsity
pattern across assemblies, making repeated contact Hessians roughly an order
of magnitude faster than the triplet path on large scenes.

.. doxygenclass:: ipc::MeshFEMHessianAssembler

Eigen Extensions
----------------

.. doxygengroup:: eigen_ext
    :content-only: