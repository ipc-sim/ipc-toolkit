if(TARGET rational::rational)
    return()
endif()

message(STATUS "Third-party: creating target 'rational::rational'")

include(CPM)
CPMAddPackage("gh:zfergus/rational-cpp#1e91438315389db5987732f464c36614b31bcadf")