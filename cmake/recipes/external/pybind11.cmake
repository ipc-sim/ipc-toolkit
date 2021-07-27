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

# Pybind11 still uses the deprecated FindPythonInterp. So let's call CMake's
# new FindPython module and set PYTHON_EXECUTABLE for Pybind11 to pick up.
# This works well with conda environments.
find_package(Python REQUIRED COMPONENTS Interpreter Development)
set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})

FetchContent_MakeAvailable(pybind11)
