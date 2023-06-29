if(TARGET GMP::GMP)
    return()
endif()

message(STATUS "Third-party: creating targets 'GMP::GMP'")

# We do not have a build recipe for this, so find it as a system installed library.
find_package(GMP REQUIRED)

if(NOT ${GMP_FOUND})
  MESSAGE(FATAL_ERROR "Unable to find GMP")
endif()

add_library(GMP_GMP INTERFACE)
add_library(GMP::GMP ALIAS GMP_GMP)

target_include_directories(GMP_GMP INTERFACE ${GMP_INCLUDE_DIRS})
target_link_libraries(GMP_GMP INTERFACE ${GMP_LIBRARIES})