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
include_guard(GLOBAL)

# options
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # When building in parallel, MSVC sometimes fails with the following error:
  # > fatal error C1090: PDB API call failed, error code '23'
  # To avoid this problem, we force PDB write to be synchronous with /FS.
  # https://developercommunity.visualstudio.com/content/problem/48897/c1090-pdb-api-call-failed-error-code-23.html
  add_compile_options(/FS)
else()
  include(ipc_toolkit_filter_flags)
  set(IPC_TOOLKIT_GLOBAL_FLAGS
      -fdiagnostics-color=always # GCC
      -fcolor-diagnostics # Clang
  )
  ipc_toolkit_filter_flags(IPC_TOOLKIT_GLOBAL_FLAGS)
  message(STATUS "Adding global flags: ${IPC_TOOLKIT_GLOBAL_FLAGS}")
  add_compile_options(${IPC_TOOLKIT_GLOBAL_FLAGS})
endif()
