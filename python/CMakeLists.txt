###############################################################################
# Create Python module for IPC Toolkit
include(pybind11)

# pybind11_add_module() adds -flto unless CMAKE_INTERPROCEDURAL_OPTIMIZATION is
# defined, and CMake passes that link flag on to nvcc's device-link step too.
# With GCC, the device-link stub then defines fatbinData in top-level asm that
# LTO cannot see, so the module fails to import with "undefined symbol:
# fatbinData". We turn LTO off for the module in CUDA builds. The tradeoff is
# small: only the binding code loses LTO, since ipc_toolkit is built without it.
# This is a normal (not cache) variable scoped to this directory, so it still
# respects a -DCMAKE_INTERPROCEDURAL_OPTIMIZATION the user passes.
if(IPC_TOOLKIT_WITH_CUDA AND NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
endif()

pybind11_add_module(ipctk)

# The compiled module is ipctk._ipctk, wrapped by the pure-Python ipctk package
# in python/ipctk/. We keep the target named ipctk so the rest of the build is
# unchanged; only the file name differs.
set_target_properties(ipctk PROPERTIES OUTPUT_NAME _ipctk)

# scikit-build-core assembles the wheel's ipctk/ folder itself: it installs
# _ipctk there (below) and copies the Python sources (wheel.packages in
# pyproject.toml). For a plain CMake build, we assemble the same package layout
# in <build>/python/ipctk/ instead, so adding <build>/python to PYTHONPATH gives
# a working `import ipctk`. We symlink the Python sources rather than copying
# them, so edits show up without a rebuild.
if(NOT SKBUILD)
  set(IPCTK_PACKAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/ipctk")
  # The generator expression stops multi-config generators from appending a
  # per-config subfolder, which would hide _ipctk from the package.
  set_target_properties(ipctk PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "$<1:${IPCTK_PACKAGE_DIR}>")

  file(MAKE_DIRECTORY "${IPCTK_PACKAGE_DIR}")
  # Linking the top-level entries is enough because subpackages link as whole
  # folders. New top-level files need a re-configure to appear.
  file(GLOB IPCTK_PYTHON_ENTRIES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/ipctk"
    "${CMAKE_CURRENT_SOURCE_DIR}/ipctk/*")
  list(FILTER IPCTK_PYTHON_ENTRIES EXCLUDE REGEX "__pycache__")
  foreach(entry IN LISTS IPCTK_PYTHON_ENTRIES)
    file(CREATE_LINK
      "${CMAKE_CURRENT_SOURCE_DIR}/ipctk/${entry}" "${IPCTK_PACKAGE_DIR}/${entry}"
      SYMBOLIC COPY_ON_ERROR)
  endforeach()
endif()

add_subdirectory(src) # Add sources to ipctk
target_include_directories(ipctk PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(ipctk PRIVATE ipc::toolkit)

# Demo simulator bindings (optional; exposed as the ipctk.demo submodule)
if(IPC_TOOLKIT_BUILD_DEMO)
  target_sources(ipctk PRIVATE "${PROJECT_SOURCE_DIR}/demo/python/simulator.cpp")
  target_link_libraries(ipctk PRIVATE ipc::toolkit::demo)
  target_compile_definitions(ipctk PRIVATE IPCTK_WITH_DEMO)

  # nlohmann::json <-> Python dict conversion for the solver parameters
  include(pybind11_json)
  target_link_libraries(ipctk PRIVATE pybind11::json)
endif()

# libigl (binding igl::edges)
include(libigl)
target_link_libraries(ipctk PRIVATE igl::core igl::predicates)

# Tight Inclusion CCD (binding tight_inclusion as ipctk.tight_inclusion)
include(tight_inclusion)
target_link_libraries(ipctk PRIVATE tight_inclusion::tight_inclusion)

# TBB (binding tbb::global_control)
include(onetbb)
target_link_libraries(ipctk PRIVATE TBB::tbb)

if (FILIB_BUILD_SHARED_LIB AND WIN32)
  # Copy DLLs to the output directory
  add_custom_command(
    TARGET ipctk POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:ipctk> $<TARGET_FILE_DIR:ipctk>
    COMMAND_EXPAND_LISTS
  )
endif()

# scikit-build-core builds the wheel from what `cmake --install` puts in the
# "python" component (see pyproject.toml), so anything the module needs at run
# time must be installed here. We place it inside the ipctk package, next to
# the Python sources that import it as ipctk._ipctk.
if(SKBUILD)
  install(TARGETS ipctk LIBRARY DESTINATION ipctk COMPONENT python)
  if(FILIB_BUILD_SHARED_LIB AND WIN32)
    # The DLLs copied above must ship next to the module too.
    install(FILES $<TARGET_RUNTIME_DLLS:ipctk> DESTINATION ipctk COMPONENT python)
  endif()
endif()

if (IPC_TOOLKIT_WITH_PROFILER)
  include(pybind11_json)
  target_link_libraries(ipctk PRIVATE pybind11::json)
endif()

# Extra warnings
# target_link_libraries(ipctk PRIVATE IPCToolkit::warnings)

# Folder name for IDE
set_target_properties(ipctk PROPERTIES FOLDER "SRC")
get_target_property(IPCTK_SOURCES ipctk SOURCES)
# The demo bindings live outside python/ (demo/python/); group them separately
list(FILTER IPCTK_SOURCES EXCLUDE REGEX "/demo/python/")
source_group(TREE "${PROJECT_SOURCE_DIR}/python" FILES ${IPCTK_SOURCES})
if(IPC_TOOLKIT_BUILD_DEMO)
  source_group("Source Files\\demo" FILES "${PROJECT_SOURCE_DIR}/demo/python/simulator.cpp")
endif()