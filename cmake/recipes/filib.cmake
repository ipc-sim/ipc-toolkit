# filib (https://github.com/zfergus/filib.git)
# License: LGPL
if(TARGET filib::filib)
  return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

include(CPM)
CPMAddPackage("gh:zfergus/filib#1cd377a7c833a68dc47217829e333f4886c5c46d")