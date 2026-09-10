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
