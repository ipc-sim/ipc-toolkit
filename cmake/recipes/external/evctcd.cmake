# Etienne Vouga's CTCD Library

if(TARGET EVCTCD)
    return()
endif()

message(STATUS "Third-party: creating target 'EVCTCD'")

include(FetchContent)
FetchContent_Declare(
    EVCTCD
    GIT_REPOSITORY https://github.com/evouga/collisiondetection.git
    GIT_TAG e5fe5c9767207df5047e375fb20180a665ae186f
    GIT_SHALLOW FALSE
)

FetchContent_GetProperties(EVCTCD)
if(NOT EVCTCD_POPULATED)
    FetchContent_Populate(EVCTCD)
endif()

file(GLOB EVCTCD_FILES "${EVCTCD_SOURCE_DIR}/src/*.cpp")
add_library(EVCTCD ${EVCTCD_FILES})
target_include_directories(EVCTCD PUBLIC "${EVCTCD_SOURCE_DIR}/include")

include(eigen)
target_link_libraries(EVCTCD PUBLIC Eigen3::Eigen)

# Turn off floating point contraction for CCD robustness
target_compile_options(EVCTCD PRIVATE "-ffp-contract=off")
