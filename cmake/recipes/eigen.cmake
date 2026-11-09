# Eigen (https://gitlab.com/libeigen/eigen, mirrored at https://github.com/eigen-mirror/eigen)
# License: MPL 2.0
if(TARGET Eigen3::Eigen)
    return()
endif()

option(EIGEN_WITH_MKL "Use Eigen with MKL" OFF)
option(EIGEN_DONT_VECTORIZE "Disable Eigen vectorization" OFF)

message(STATUS "Third-party: creating target 'Eigen3::Eigen'")

include(CPM)

set(EIGEN_GIT_TAG 5.0.1)
set(EIGEN_GITLAB_URL "https://gitlab.com/libeigen/eigen.git")
set(EIGEN_GITHUB_URL "https://github.com/eigen-mirror/eigen.git")

# Cloning from GitLab fails often enough in CI that we fall back to the GitHub
# mirror, which carries the same tags and commits. A failed clone inside
# CPMAddPackage is a fatal error we can't recover from, so instead we probe
# GitLab up front with a cheap `git ls-remote` and pick the source before CPM
# runs. We skip the probe whenever no download will happen: Eigen is already in
# the CPM cache or this build's _deps, or the user pointed CPM at a local copy.
if(CPM_SOURCE_CACHE)
    set(EIGEN_DOWNLOAD_DIR "${CPM_SOURCE_CACHE}/eigen/${EIGEN_GIT_TAG}")
else()
    set(EIGEN_DOWNLOAD_DIR "${FETCHCONTENT_BASE_DIR}/eigen-src")
endif()
set(EIGEN_GIT_REPOSITORY "${EIGEN_GITLAB_URL}")
# CPM only looks for git once it starts downloading, so we look for it here.
find_package(Git QUIET)
if(GIT_FOUND AND NOT (EXISTS "${EIGEN_DOWNLOAD_DIR}"
        OR DEFINED CPM_eigen_SOURCE OR FETCHCONTENT_FULLY_DISCONNECTED))
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" ls-remote --exit-code
            "${EIGEN_GITLAB_URL}" "refs/tags/${EIGEN_GIT_TAG}"
        RESULT_VARIABLE EIGEN_GITLAB_STATUS
        OUTPUT_QUIET ERROR_QUIET
        TIMEOUT 30
    )
    if(NOT EIGEN_GITLAB_STATUS EQUAL 0)
        message(WARNING "Unable to reach ${EIGEN_GITLAB_URL}; "
            "falling back to the GitHub mirror ${EIGEN_GITHUB_URL}.")
        set(EIGEN_GIT_REPOSITORY "${EIGEN_GITHUB_URL}")
    endif()
endif()

CPMAddPackage(
    NAME eigen
    GIT_REPOSITORY "${EIGEN_GIT_REPOSITORY}"
    GIT_TAG ${EIGEN_GIT_TAG}
    # Both sources hold identical commits, so we key the cache on the tag alone.
    # Otherwise CPM hashes the repository URL into the cache path and a
    # fallback would re-clone Eigen even when the GitLab copy is already cached.
    CUSTOM_CACHE_KEY ${EIGEN_GIT_TAG}
    DOWNLOAD_ONLY YES
)

add_library(Eigen3_Eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS Eigen3_Eigen)

include(GNUInstallDirs)
target_include_directories(Eigen3_Eigen SYSTEM INTERFACE
    $<BUILD_INTERFACE:${eigen_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties(Eigen3_Eigen PROPERTIES
  INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}"
)

if(EIGEN_DONT_VECTORIZE)
    target_compile_definitions(Eigen3_Eigen INTERFACE EIGEN_DONT_VECTORIZE=1)
endif()

if(EIGEN_WITH_MKL)
    # TODO: Checks that, on 64bits systems, `mkl::mkl` is using the LP64 interface
    # (by looking at the compile definition of the target)
    include(mkl)
    target_link_libraries(Eigen3_Eigen INTERFACE mkl::mkl)
    target_compile_definitions(Eigen3_Eigen INTERFACE
        EIGEN_USE_MKL_ALL
        EIGEN_USE_LAPACKE_STRICT
    )
elseif(APPLE)
    find_package(BLAS REQUIRED)
    find_library(LAPACKE lapacke PATHS
        "/opt/local/lib/lapack"
        "/opt/homebrew/opt/lapack/lib"
    )
    if (NOT LAPACKE)
        # BLAS should be available on macOS, but LAPACKE might not be
        message(WARNING "LAPACKE library not found (required for EIGEN_USE_LAPACKE on macOS)! "
            "Perhaps you need to install it (e.g., brew install lapack). "
            "Eigen will be built without LAPACKE support.")
    else()
        message(STATUS "Found BLAS and LAPACKE. Enabling Eigen LAPACKE support.")
        target_link_libraries(Eigen3_Eigen INTERFACE
            ${BLAS_LIBRARIES} ${LAPACKE}
        )
        target_compile_definitions(Eigen3_Eigen INTERFACE
            EIGEN_USE_BLAS
            EIGEN_USE_LAPACKE_STRICT
        )
    endif()
endif()

# On Windows, enable natvis files to improve debugging experience
if(WIN32 AND eigen_SOURCE_DIR)
    target_sources(Eigen3_Eigen INTERFACE $<BUILD_INTERFACE:${eigen_SOURCE_DIR}/debug/msvc/eigen.natvis>)
endif()

# Install rules
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME eigen)
set_target_properties(Eigen3_Eigen PROPERTIES EXPORT_NAME Eigen)
install(DIRECTORY ${eigen_SOURCE_DIR} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS Eigen3_Eigen EXPORT Eigen_Targets)
install(EXPORT Eigen_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/eigen NAMESPACE Eigen3::)