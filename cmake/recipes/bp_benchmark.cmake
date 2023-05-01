if(TARGET bp_benchmark)
    return()
endif()

message(STATUS "Third-party: creating target 'bp_benchmark'")

include(FetchContent)
FetchContent_Declare(
    bp_benchmark
    GIT_REPOSITORY https://github.com/dbelgrod/broad-phase-benchmark.git
    GIT_TAG main
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(bp_benchmark)

add_library(bp_benchmark ALIAS BroadPhaseBenchmark)