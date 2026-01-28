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

Is there a way to ignore select collisions?
-------------------------------------------

Yes, it is possible to ignore select collisions.

The functionality for doing so is through the :cpp:`BroadPhase::can_vertices_collide`.
This function takes two vertex IDs and returns a true if the vertices can collide otherwise false.

This is used to determine if any geometry connected to the verties can collide. E.g., when checking if vertex ``vi`` can collide with triangle ``f = (vj, vk, vl)``, the code checks:

.. code-block::

    can_face_vertex_collide(f, vi) := can_vertices_collide(vj, vi) && can_vertices_collide(vk, vi) && can_vertices_collide(vl, vi)

This is a little limited since it will ignore the one-ring around a vertex instead of a single face-vertex pair, but hopefully that can get you started.

To get something more customized, you can try to modify the BroadPhase class, which has these functions hard-coded:

.. code-block:: c++

    virtual bool can_edge_vertex_collide(size_t ei, size_t vi) const;
    virtual bool can_edges_collide(size_t eai, size_t ebi) const;
    virtual bool can_face_vertex_collide(size_t fi, size_t vi) const;
    virtual bool can_edge_face_collide(size_t ei, size_t fi) const;
    virtual bool can_faces_collide(size_t fai, size_t fbi) const;

You can modify these with function pointers or override them to have the specific implementation you are interested in.

.. note::

    If you are building collisions through the ``Candidates`` class, the ``Candidates::build`` function sets the ``BroadPhase::can_vertices_collide`` using the ``CollisionMesh::can_collide`` function pointer. This ``CollisionMesh::can_collide`` function uses the same interface as the ``BroadPhase::can_vertices_collide`` above.

.. warning::
    This method is not recommended for Python since calling a Python lambda function from the C++ side is too slow to use. Instead there are ``SparseCanCollide`` and ``VertexPatchesCanCollide`` classes in Python to help do this efficiently.

My question is not answered here. What should I do?
---------------------------------------------------

Please open an Q&A discussion post on `GitHub <https://github.com/ipc-sim/ipc-toolkit/discussions>`_ and we will do our best to help you.