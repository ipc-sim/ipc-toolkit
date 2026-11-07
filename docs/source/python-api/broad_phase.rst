Broad Phase
===========

Broad Phase
-----------

.. autoclass:: ipctk.BroadPhase

    .. autoclasstoc::

Brute Force
-----------

.. autoclass:: ipctk.BruteForce

    .. autoclasstoc::

Hash Grid
---------

.. autoclass:: ipctk.HashGrid

    .. autoclasstoc::

Spatial Hash
------------

.. autoclass:: ipctk.SpatialHash

    .. autoclasstoc::

LBVH
----

.. autoclass:: ipctk.LBVH

    .. autoclasstoc::

Sweep and Prune
---------------

.. autoclass:: ipctk.SweepAndPrune

    .. autoclasstoc::

Sweep and Tiniest Queue
-----------------------

``ipctk.SweepAndTiniestQueue`` is available only when ``ipctk`` is built with
CUDA (``IPC_TOOLKIT_WITH_CUDA``), which the documentation build is not.

.. .. autoclass:: ipctk.SweepAndTiniestQueue
..
..     .. autoclasstoc::

LBVH (CUDA)
-----------

``ipctk.cuda.LBVH`` is the GPU counterpart of ``ipctk.LBVH`` (C++:
``ipc::cuda::LBVH``). The ``ipctk.cuda`` submodule mirrors the C++ ``ipc::cuda``
namespace; its classes exist only when ``ipctk`` is built with CUDA
(``IPC_TOOLKIT_WITH_CUDA``), which the documentation build is not. The CUDA
wheels also import on machines without a GPU, where using ``ipctk.cuda.LBVH``
raises; check ``ipctk.cuda.is_available()`` first.

.. .. autoclass:: ipctk.cuda.LBVH
..
..     .. autoclasstoc::

AABB
----

.. autoclass:: ipctk.AABB

    .. autoclasstoc::