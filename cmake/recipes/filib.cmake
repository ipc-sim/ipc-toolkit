# filib (https://github.com/zfergus/filib.git)
# License: LGPL
if(TARGET filib::filib)
  return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

include(CPM)
CPMAddPackage("gh:zfergus/filib#63acb1f78dcdf7667ca2d214f7a3fa3505fdf339")