# autodiff (https://github.com/autodiff/autodiff.git)
# BSD license
if(TARGET autodiff::autodiff)
    return()
endif()

message(STATUS "Third-party: creating target 'autodiff::autodiff'")

option(AUTODIFF_BUILD_TESTS "Enable the compilation of the test files." OFF)
option(AUTODIFF_BUILD_PYTHON "Enable the compilation of the python bindings." OFF)
option(AUTODIFF_BUILD_EXAMPLES "Enable the compilation of the example files." OFF)
option(AUTODIFF_BUILD_DOCS "Enable the build of the documentation and website." OFF)

include(CPM)
CPMAddPackage("gh:autodiff/autodiff@1.0.3")
target_link_libraries(ipc_toolkit PRIVATE autodiff::autodiff)
