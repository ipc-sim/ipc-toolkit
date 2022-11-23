if(TARGET filib::filib)
    return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

include(FetchContent)
FetchContent_Declare(
    filib
    GIT_REPOSITORY https://github.com/txstc55/filib.git
    GIT_TAG 31c801fc45b545a7809911b193498efdb4bd930d
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(filib)

add_library(filib::filib ALIAS filib)
