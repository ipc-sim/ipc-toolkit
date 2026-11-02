# CCD Query IO (https://github.com/Continuous-Collision-Detection/CCD-Query-IO)
# License: MIT
if(TARGET ccd_io::ccd_io)
    return()
endif()

message(STATUS "Third-party: creating target 'ccd_io::ccd_io'")

include(ipc_toolkit_tests_data)

include(CPM)
CPMAddPackage(
    URI "gh:Continuous-Collision-Detection/CCD-Query-IO#fbead5891e521868d5d449eb656cbcf70c13d5d3"
    OPTIONS "CCD_IO_DOWNLOAD_SAMPLE_QUERIES ON"
            "CCD_IO_SAMPLE_QUERIES_DIR ${IPC_TOOLKIT_TESTS_DATA_DIR}/ccd-queries/"
)