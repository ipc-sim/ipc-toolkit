# CCD Query IO (https://github.com/Continuous-Collision-Detection/CCD-Query-IO)
# License: MIT
if(TARGET ccd_io::ccd_io)
    return()
endif()

message(STATUS "Third-party: creating target 'ccd_io::ccd_io'")

include(ipc_toolkit_tests_data)

include(CPM)
CPMAddPackage(
    URI "gh:Continuous-Collision-Detection/CCD-Query-IO#36f6093af81a65acc27d9f05ad32d6b5729e8d15"
    OPTIONS "CCD_IO_DOWNLOAD_SAMPLE_QUERIES ON"
            "CCD_IO_SAMPLE_QUERIES_DIR ${IPC_TOOLKIT_TESTS_DATA_DIR}/ccd-queries/"
)