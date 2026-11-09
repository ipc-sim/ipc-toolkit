Utils
=====

Logger
------

.. doxygenfunction:: ipc::logger
.. doxygenfunction:: ipc::set_logger

CUDA Availability
-----------------

.. doxygenfunction:: ipc::cuda::is_available

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

Eigen ABI Check
---------------

Eigen picks its alignment, and so the size and layout of fixed-size types and
the way it allocates memory, from the flags each translation unit is compiled
with. The IPC Toolkit defines ``EIGEN_DONT_VECTORIZE`` when SIMD is enabled, so
code that passes Eigen objects to it must be compiled the same way. Linking the
``ipc::toolkit`` CMake target does this for you. Use these functions to check it
at runtime, for example in a test.

.. doxygenstruct:: ipc::EigenAbi
    :members:
.. doxygenfunction:: ipc::eigen_abi
.. doxygenfunction:: ipc::library_eigen_abi
.. doxygenfunction:: ipc::eigen_abi_matches_library

Eigen Extensions
----------------

.. doxygengroup:: eigen_ext
    :content-only: