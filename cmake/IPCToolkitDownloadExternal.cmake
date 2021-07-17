include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
  set(IPC_TOOLKIT_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
  set(IPC_TOOLKIT_EXTRA_OPTIONS "")
endif()

function(ipc_toolkit_download_project name)
  download_project(
    PROJ         ${name}
    SOURCE_DIR   ${IPC_TOOLKIT_EXTERNAL}/${name}
    DOWNLOAD_DIR ${IPC_TOOLKIT_EXTERNAL}/.cache/${name}
    QUIET
    ${IPC_TOOLKIT_EXTRA_OPTIONS}
    ${ARGN}
  )
endfunction()

################################################################################

function(ipc_toolkit_download_catch2)
  ipc_toolkit_download_project(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.12.2
  )
endfunction()

function(ipc_toolkit_download_libigl)
  ipc_toolkit_download_project(libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG        v2.3.0
  )
endfunction()

function(ipc_toolkit_download_finite_diff)
  ipc_toolkit_download_project(finite-diff
    GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
    GIT_TAG        3cc4537130708fa08ce2a4cd0d77b2621227ccfe
  )
endfunction()

function(ipc_toolkit_download_tbb)
  ipc_toolkit_download_project(tbb
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        141b0e310e1fb552bdca887542c9c1a8544d6503
  )
endfunction()

# Etienne Vouga's CTCD Library
function(ipc_toolkit_download_evctcd)
  ipc_toolkit_download_project(EVCTCD
    GIT_REPOSITORY https://github.com/evouga/collisiondetection.git
    GIT_TAG        e5fe5c9767207df5047e375fb20180a665ae186f
  )
endfunction()

# Provably conservative CCD of [Wang and Ferguson et al. 2020]
function(ipc_toolkit_download_tight_inclusion)
  ipc_toolkit_download_project(Tight-Inclusion
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Tight-Inclusion.git
    GIT_TAG        90450c3690c06e116056f225c3b7c0ddef77b54f
  )
endfunction()

function(ipc_toolkit_download_spdlog)
  ipc_toolkit_download_project(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.8.0
  )
endfunction()

function(ipc_toolkit_download_fmt)
  ipc_toolkit_download_project(fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG        8.0.1
  )
endfunction()

function(ipc_toolkit_download_robin_map)
  ipc_toolkit_download_project(robin-map
    GIT_REPOSITORY https://github.com/Tessil/robin-map.git
    GIT_TAG        v0.6.3
  )
endfunction()
