if(TARGET STQ::CPU)
    return()
endif()

message(STATUS "Third-party: creating target 'STQ::CPU'")

option(STQ_WITH_CPU  "Enable CPU Implementation"   ON)
option(STQ_WITH_CUDA "Enable CUDA Implementation" OFF) # get this through GPU_CCD

    include(FetchContent)
    FetchContent_Declare(
        sweep_and_tiniest_queue
        GIT_REPOSITORY https://github.com/dbelgrod/broadphase-gpu.git
        GIT_TAG 8ea831260de6c0100b597f7568434cd4e6eac580
        GIT_SHALLOW FALSE
    )
    FetchContent_MakeAvailable(sweep_and_tiniest_queue)