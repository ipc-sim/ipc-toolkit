# Etienne Vouga's CCD Library (https://github.com/evouga/collisiondetection.git)
# License: ???
if(TARGET evouga::ccd)
    return()
endif()

message(STATUS "Third-party: creating target 'evouga::ccd'")

include(CPM)
CPMAddPackage(
    NAME evccd
    GITHUB_REPOSITORY evouga/collisiondetection
    GIT_TAG e5fe5c9767207df5047e375fb20180a665ae186f
    DOWNLOAD_ONLY YES
)

# file(GLOB EVOUGA_CCD_SOURCE_FILES "${evccd_SOURCE_DIR}/src/*.cpp")
add_library(evouga_ccd
    "${evccd_SOURCE_DIR}/src/CTCD.cpp"
)
add_library(evouga::ccd ALIAS evouga_ccd)

target_include_directories(evouga_ccd PUBLIC "${evccd_SOURCE_DIR}/include")

include(eigen)
target_link_libraries(evouga_ccd PUBLIC Eigen3::Eigen)

# Turn off floating point contraction for CCD robustness
target_compile_options(evouga_ccd PRIVATE "-ffp-contract=off")
