# ipc-toolkit-test-data (https://github.com/ipc-sim/ipc-toolkit-test-data)
# License: MIT

if(TARGET ipc_toolkit_test_data_download)
    return()
endif()

include(ExternalProject)
include(FetchContent)

set(IPC_TOOLKIT_TESTS_DATA_DIR "${PROJECT_SOURCE_DIR}/tests/data/" CACHE PATH "Where should we download the tests data?")
option(IPC_TOOLKIT_USE_EXISTING_TESTS_DATA_DIR "Use an existing test data directory instead of downloading it" OFF)

if(IPC_TOOLKIT_USE_EXISTING_TESTS_DATA_DIR)
    ExternalProject_Add(
        ipc_toolkit_test_data_download
        PREFIX "${FETCHCONTENT_BASE_DIR}/tests/data"
        SOURCE_DIR ${IPC_TOOLKIT_TESTS_DATA_DIR}

        # NOTE: No download step
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
    )
else()
    ExternalProject_Add(
        ipc_toolkit_test_data_download
        PREFIX "${FETCHCONTENT_BASE_DIR}/tests/data"
        SOURCE_DIR ${IPC_TOOLKIT_TESTS_DATA_DIR}

        GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit-tests-data.git
        GIT_TAG 9061bb85b30b7e410db736f6e16de234d5b0555a

        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
        LOG_DOWNLOAD ON
    )
endif()