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
  if(NOT TARGET tbb_static)
    ipc_toolkit_download_tbb()
    set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
    set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
    set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
    add_subdirectory(${IPC_TOOLKIT_EXTERNAL}/tbb EXCLUDE_FROM_ALL)
  endif()
  add_library(TBB::tbb ALIAS tbb_static)
endif()

# CCD
if(IPC_TOOLKIT_WITH_CORRECT_CCD)
  # Provably conservative CCD of [Wang and Ferguson et al. 2020]
  if(NOT TARGET TightInclusion)
    ipc_toolkit_download_tight_inclusion()
    add_subdirectory(${IPC_TOOLKIT_EXTERNAL}/Tight-Inclusion)
    add_library(TightInclusion ALIAS tight_inclusion)
  endif()
else()
  # Etienne Vouga's CTCD Library
  if(NOT TARGET EVCTCD)
    ipc_toolkit_download_evctcd()

    file(GLOB EVCTCD_FILES "${IPC_TOOLKIT_EXTERNAL}/EVCTCD/src/*.cpp")
    add_library(EVCTCD ${EVCTCD_FILES})
    target_include_directories(EVCTCD PUBLIC "${IPC_TOOLKIT_EXTERNAL}/EVCTCD/include")
    target_link_libraries(EVCTCD PUBLIC Eigen3::Eigen)

    # Turn off floating point contraction for CCD robustness
    target_compile_options(EVCTCD PUBLIC "-ffp-contract=off")
  endif()
endif()


# Logger
if(IPC_TOOLKIT_WITH_LOGGER AND NOT TARGET spdlog::spdlog)
  ipc_toolkit_download_spdlog()
  add_library(spdlog INTERFACE)
  add_library(spdlog::spdlog ALIAS spdlog)
  target_include_directories(spdlog SYSTEM INTERFACE ${IPC_TOOLKIT_EXTERNAL}/spdlog/include)
endif()
