# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(IPCToolkitDownloadExternal)

################################################################################
# Required libraries
################################################################################

# libigl
if(NOT TARGET igl::core)
    ipc_toolkit_download_libigl()
    # Import libigl targets
    list(APPEND CMAKE_MODULE_PATH "${IPC_TOOLKIT_EXTERNAL}/libigl/cmake")
    include(libigl)
endif()

# TBB
if(NOT TARGET TBB::tbb)
  ipc_toolkit_download_tbb()
  set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
  set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  add_subdirectory(${IPC_TOOLKIT_EXTERNAL}/tbb EXCLUDE_FROM_ALL)
  add_library(TBB::tbb ALIAS tbb_static)
endif()

# finite-diff
if(NOT TARGET FiniteDiff::FiniteDiff)
  ipc_toolkit_download_finite_diff()
  add_subdirectory(${IPC_TOOLKIT_EXTERNAL}/finite-diff EXCLUDE_FROM_ALL)
  add_library(FiniteDiff::FiniteDiff ALIAS FiniteDiff)
endif()

# Etienne Vouga's CTCD Library
if(NOT TARGET EVCTCD)
  ipc_toolkit_download_evctcd()
  # Set Eigen directory environment variable (needed for EVCTCD)
  set(ENV{EIGEN3_INCLUDE_DIR} "${IPC_TOOLKIT_EXTERNAL}/libigl/external/eigen/")
  add_subdirectory(${IPC_TOOLKIT_EXTERNAL}/EVCTCD)
  # These includes are PRIVATE for some reason
  target_include_directories(collisiondetection PUBLIC "${IPC_TOOLKIT_EXTERNAL}/EVCTCD/include")
  # Turn off floating point contraction for CCD robustness
  target_compile_options(collisiondetection PUBLIC "-ffp-contract=off")
  # Rename for convenience
  add_library(EVCTCD ALIAS collisiondetection)
endif()

# Rational implmentation of Brochu et al. [2012]
if(NOT TARGET RationalCCD)
  ipc_toolkit_download_rational_ccd()
  add_subdirectory(${IPC_TOOLKIT_EXTERNAL}/rational_ccd)
endif()
