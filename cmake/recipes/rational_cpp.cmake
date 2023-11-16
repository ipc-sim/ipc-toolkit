# rational-cpp (https://github.com/zfergus/rational-cpp)
# License: MIT
if(TARGET rational::rational)
    return()
endif()

message(STATUS "Third-party: creating target 'rational::rational'")

include(CPM)
CPMAddPackage("gh:zfergus/rational-cpp#687d4ea3436ada7231b8920f3cd5b02b438c21aa")