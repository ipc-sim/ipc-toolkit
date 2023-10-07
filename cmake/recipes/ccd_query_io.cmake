# CCD Query IO (https://github.com/Continuous-Collision-Detection/CCD-Query-IO)
# License: MIT
if(TARGET ccd_io::ccd_io)
    return()
endif()

message(STATUS "Third-party: creating target 'ccd_io::ccd_io'")

set(CCD_IO_DOWNLOAD_SAMPLE_QUERIES ON CACHE BOOL "Download sample CCD queries" FORCE)
set(CCD_IO_SAMPLE_QUERIES_DIR "${PROJECT_SOURCE_DIR}/tests/data/ccd-queries/" CACHE PATH "Where should we download sample queries?")

include(CPM)
CPMAddPackage("gh:Continuous-Collision-Detection/CCD-Query-IO#36f6093af81a65acc27d9f05ad32d6b5729e8d15")