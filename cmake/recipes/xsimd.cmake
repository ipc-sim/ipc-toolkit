# xsimd (https://github.com/xtensor-stack/xsimd)
# License: BSD-3-Clause
if(TARGET xsimd::xsimd)
  return()
endif()

message(STATUS "Third-party: creating target 'xsimd::xsimd'")

include(CPM)
CPMAddPackage("gh:xtensor-stack/xsimd#14.0.0")

add_library(xsimd::xsimd ALIAS xsimd)

# Folder name for IDE
set_target_properties(xsimd PROPERTIES FOLDER "ThirdParty")
