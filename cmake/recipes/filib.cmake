# filib (https://github.com/zfergus/filib.git)
# License: LGPL
if(TARGET filib::filib)
  return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

include(CPM)
CPMAddPackage("gh:zfergus/filib#e09f00eb20e8cf04ac5cfc109f69560284c4e51a")