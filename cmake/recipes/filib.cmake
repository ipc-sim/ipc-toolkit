# filib (https://github.com/zfergus/filib.git)
# License: LGPL
if(TARGET filib::filib)
  return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

include(CPM)
CPMAddPackage("gh:zfergus/filib#636f54db74bf4f658779f334ba04aeef96b8d392")