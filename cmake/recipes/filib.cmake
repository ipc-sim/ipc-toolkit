# filib (https://github.com/zfergus/filib.git)
# License: LGPL
if(TARGET filib::filib)
  return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

include(CPM)
CPMAddPackage("gh:zfergus/filib#fe1d1bfa56ef8c00fbcea2deec69a6b4f9da2207")