.. _faq:

Frequently Asked Questions
==========================

.. role:: cpp(code)
   :language: c++
.. role:: cmake(code)
   :language: cmake
.. role:: python(code)
   :language: python

How do I include IPC Toolkit in my project?
-------------------------------------------

If you are using CMake, the public include directory is added to the `ipc::toolkit` cmake target which means that any lib/bin that includes `ipc::toolkit` as a dependency also adds those include directories too.

If you are not using CMake, the include path is ``src``.

Files are included with the prefix :cpp:`#include <ipc/...>` in C++ and :python:`import ipctk` in Python.

How do I determine which edges intersect?
-----------------------------------------

We do not provide an edge-edge intersection function in 3D, but you can approximate it by computing the distance between the two edges and checking if it is less than a threshold.

Inside any :cpp:`BroadPhase` classes, the function :cpp:`detect_edge_edge_candidates` determines which edges intersect based on their bounding boxes. You can then use :cpp:`edge_edge_distance` to check if they approximatly intersect by computing the distance.

How do I build the edge matrix from the face matrix?
----------------------------------------------------

To build the edge matrix you can use :cpp:`igl::edges(faces, edges);` in C++ or :python:`ipctk.edges(faces)` in Python.

My question is not answered here. What should I do?
---------------------------------------------------

Please open an issue on `GitHub <https://github.com/ipc-sim/ipc-toolkit/issues>`_ and we will do our best to help you.