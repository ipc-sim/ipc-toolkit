if(TARGET ipc::toolkit::test_data)
    return()
endif()

include(ExternalProject)
include(FetchContent)

set(IPC_TOOLKIT_TEST_DATA_DIR "${PROJECT_SOURCE_DIR}/tests/data/" CACHE PATH "Where should we download the test data?")

ExternalProject_Add(
    ipc_toolkit_test_data_download
    PREFIX "${FETCHCONTENT_BASE_DIR}/tests/data"
    SOURCE_DIR ${IPC_TOOLKIT_TEST_DATA_DIR}

    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit-test-data.git
    GIT_TAG a4fe14a4fd689a013ae25b72e50bdccac3bf2d04

    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
)

# Create a dummy target for convenience
add_library(ipc_toolkit_test_data INTERFACE)
add_library(ipc::toolkit::test_data ALIAS ipc_toolkit_test_data)

add_dependencies(ipc_toolkit_test_data ipc_toolkit_test_data_download)

target_compile_definitions(ipc_toolkit_test_data INTERFACE IPC_TOOLKIT_TEST_DATA_DIR=\"${IPC_TOOLKIT_TEST_DATA_DIR}\")