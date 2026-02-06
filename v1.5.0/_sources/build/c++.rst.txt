Building Library
================

.. role:: cpp(code)
   :language: c++
.. role:: cmake(code)
   :language: cmake

IPC Toolkit uses CMake to configure its build system. It is designed to be used as a submodule in your project, and as such does not have a proper install target. Instead, you can include it in your project using CMake's ``FetchContent`` module.

Adding IPC Toolkit to your CMake
--------------------------------

The easiest way to add the toolkit to an existing CMake project is to download it through CMake.
CMake provides functionality for doing this called `FetchContent <https://cmake.org/cmake/help/latest/module/FetchContent.html>`_ (requires CMake ≥ 3.14).
We use a very similar process to download all external dependencies (using `CPM.cmake <https://github.com/cpm-cmake/CPM.cmake>`_).

For example,

.. code-block:: cmake

   include(FetchContent)
   FetchContent_Declare(
       ipc_toolkit
       GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
       GIT_TAG ${IPC_TOOLKIT_GIT_TAG}
   )
   FetchContent_MakeAvailable(ipc_toolkit)

where ``IPC_TOOLKIT_GIT_TAG`` is set to the version of the toolkit you want to use. This will download and add the toolkit to CMake. The toolkit can then be linked against using

.. code-block:: cmake

   # Link against the IPC Toolkit
   target_link_libraries(${PROJECT_NAME} PUBLIC ipc::toolkit)

where ``PROJECT_NAME`` is the name of your library/binary.

.. tip::
   If your ``IPC_TOOLKIT_GIT_TAG`` is a tag (e.g. ``v1.3.1``), then you can use the ``FetchContent_Declare`` argument ``GIT_SHALLOW TRUE`` to download only a single commit. Otherwise, you should use the default ``GIT_SHALLOW FALSE``.

Building the Library
--------------------

You can build the IPC Toolkit using CMake as you would any other CMake project. The following is a minimal example of how to build the toolkit:

.. code-block:: bash

   mkdir build
   cd build
   cmake -DIPC_TOOLKIT_BUILD_TESTS=ON -DIPC_TOOLKIT_BUILD_PYTHON=ON ..
   make -j$(nproc)

This will build the IPC Toolkit and all of its dependencies. The ``IPC_TOOLKIT_BUILD_TESTS`` option enables building the unit tests, and the ``IPC_TOOLKIT_BUILD_PYTHON`` option enables building the Python bindings.

CUDA support is disabled by default. Enable it by setting CMake option ``IPC_TOOLKIT_WITH_CUDA`` to ON.

.. warning::
   Installing the IPC Toolkit using the ``make install`` has not been tested and is not recommended. The IPC Toolkit is designed to be used as a submodule in your project, and as such does not have a proper install target.

Dependencies
------------

**All required dependencies are downloaded through CMake** depending on the build options, and are built automatically when you build the IPC Toolkit. You do not need to install them separately.

A full list of dependencies can be found on the `dependencies page <https://ipctk.xyz/dependencies.html>`_.
