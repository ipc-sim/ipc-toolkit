################################################################################
# IPC Toolkit Demo Simulator
#
# Application-side code intentionally kept out of the core library: the
# incremental potential presented to the (polysolve) nonlinear solver and the
# time-stepping loop. The core library only provides the modular components
# (energies, time integrators, potentials, and CCD); this target is the living
# proof that they are easy to assemble into a simulator.
################################################################################

set(SOURCES
  src/kinematics.cpp
  src/kinematics.hpp
  src/simulator.cpp
  src/simulator.hpp
)

add_library(ipc_toolkit_demo ${SOURCES})
add_library(ipc::toolkit::demo ALIAS ipc_toolkit_demo)

target_include_directories(ipc_toolkit_demo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

# The demo depends on the core library; the core library must never depend on
# the demo (or on any solver).
target_link_libraries(ipc_toolkit_demo PUBLIC ipc::toolkit)

# polysolve (nonlinear solver with PolyFEM-grade Newton + line search).
# PUBLIC: the Simulator IS a polysolve::nonlinear::Problem.
include(polysolve)
target_link_libraries(ipc_toolkit_demo PUBLIC polysolve::polysolve)

target_compile_features(ipc_toolkit_demo PUBLIC cxx_std_17)

# Folder name for IDE
set_target_properties(ipc_toolkit_demo PROPERTIES FOLDER "DEMO")
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${SOURCES})

# NOTE: The demo's tests live in the unified ipc_toolkit_tests binary
# (tests/src/tests/demo/) and its Python bindings in the core ipctk module as
# the ipctk.demo submodule (demo/python/simulator.cpp); both are added by
# their respective parent targets when IPC_TOOLKIT_BUILD_DEMO is ON.
