#
# Copyright 2020 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#
if(TARGET pybind11::pybind11)
    return()
endif()

message(STATUS "Third-party: creating target 'pybind11::pybind11'")

include(FetchContent)
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.7.0
    GIT_SHALLOW TRUE
)

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

# Pybind11 still uses the deprecated FindPythonInterp. So let's call CMake's
# new FindPython module and set PYTHON_EXECUTABLE for Pybind11 to pick up.
# This works well with conda environments.
find_package(Python REQUIRED COMPONENTS Interpreter Development)
set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})

FetchContent_MakeAvailable(pybind11)
