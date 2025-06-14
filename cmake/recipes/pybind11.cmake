# pybind11 (https://github.com/pybind/pybind11)
# License: BSD-style
if(TARGET pybind11::pybind11)
    return()
endif()

message(STATUS "Third-party: creating target 'pybind11::pybind11'")

if (POLICY CMP0094)  # https://cmake.org/cmake/help/latest/policy/CMP0094.html
    cmake_policy(SET CMP0094 NEW)  # FindPython should return the first matching Python
endif ()

# needed on GitHub Actions CI: actions/setup-python does not touch registry/frameworks on Windows/macOS
# this mirrors PythonInterp behavior which did not consult registry/frameworks first
if (NOT DEFINED Python_FIND_REGISTRY)
    set(Python_FIND_REGISTRY "LAST")
endif ()
if (NOT DEFINED Python_FIND_FRAMEWORK)
    set(Python_FIND_FRAMEWORK "LAST")
endif ()

include(CPM)
CPMAddPackage("gh:pybind/pybind11@2.13.1")
