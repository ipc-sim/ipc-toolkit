#
# Copyright 2021 Adobe. All rights reserved.
# This file is licensed to you under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License. You may obtain a copy
# of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under
# the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR REPRESENTATIONS
# OF ANY KIND, either express or implied. See the License for the specific language
# governing permissions and limitations under the License.
#
function(ipc_toolkit_filter_flags flags)
  include(CheckCXXCompilerFlag)
  set(output_flags)
  foreach(FLAG IN ITEMS ${${flags}})
    string(REPLACE "=" "-" FLAG_VAR "${FLAG}")
    if(NOT DEFINED IS_SUPPORTED_${FLAG_VAR})
      check_cxx_compiler_flag("${FLAG}" IS_SUPPORTED_${FLAG_VAR})
    endif()
    if(IS_SUPPORTED_${FLAG_VAR})
      list(APPEND output_flags $<$<COMPILE_LANGUAGE:CXX>:${FLAG}>)
    endif()
  endforeach()
  set(${flags} ${output_flags} PARENT_SCOPE)
endfunction()

# Check that the C++ compiler accepts all of `flags` together, and set `result`
# to TRUE or FALSE in the caller's scope. Both arguments are variable names.
#
# Unlike ipc_toolkit_filter_flags(), which tests each flag on its own and drops
# the ones it rejects, this tests the whole list at once and reports a single
# answer. We use it when the flags were chosen by the user, where silently
# dropping a flag would hide a mistake.
#
# We only compile a trivial program and never run it, so flags that target a
# different CPU than the one doing the build still pass. The answer is cached
# under a name built from the flags, so changing the flags re-runs the check
# instead of reusing an old result.
function(ipc_toolkit_check_cxx_flags flags result)
  include(CheckCXXSourceCompiles)
  string(MAKE_C_IDENTIFIER "${${flags}}" flags_id)
  # CMAKE_REQUIRED_FLAGS wants one space-separated string, not a CMake list. We
  # set it inside this function so it does not leak into the caller's scope.
  string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${${flags}}")
  check_cxx_source_compiles(
    "int main() { return 0; }" IPC_TOOLKIT_CXX_FLAGS_OK${flags_id})
  set(${result} ${IPC_TOOLKIT_CXX_FLAGS_OK${flags_id}} PARENT_SCOPE)
endfunction()

# The nvcc counterpart of ipc_toolkit_filter_flags(): keep the flags of `flags`
# that nvcc's host compiler accepts, checked by actually compiling with nvcc and
# `-Xcompiler=<flag>`, and wrap each so it applies only to CUDA sources compiled
# by nvcc. `-Xcompiler` is required: nvcc parses some host flags itself with a
# different meaning (`-Werror` takes nvcc's own diagnostic names, `-march=...`
# is read as an input file), so a bare host flag on the nvcc command line is
# unsafe. Requires the CUDA language to be enabled.
function(ipc_toolkit_filter_nvcc_flags flags)
  include(CheckCompilerFlag)
  set(output_flags)
  foreach(FLAG IN ITEMS ${${flags}})
    string(REPLACE "=" "-" FLAG_VAR "${FLAG}")
    if(NOT DEFINED IS_SUPPORTED_NVCC_HOST_${FLAG_VAR})
      check_compiler_flag(CUDA "-Xcompiler=${FLAG}" IS_SUPPORTED_NVCC_HOST_${FLAG_VAR})
    endif()
    if(IS_SUPPORTED_NVCC_HOST_${FLAG_VAR})
      list(APPEND output_flags
        "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CUDA_COMPILER_ID:NVIDIA>>:-Xcompiler=${FLAG}>")
    endif()
  endforeach()
  set(${flags} ${output_flags} PARENT_SCOPE)
endfunction()
