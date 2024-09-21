# filib (https://github.com/zfergus/filib.git)
# License: LGPL-2.1
if(TARGET filib::filib)
  return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

# This has to be set to ON to avoid licensing IPC Toolkit under LGPL
set(FILIB_BUILD_SHARED_LIB ON CACHE BOOL "Build shared library" FORCE)

include(CPM)
CPMAddPackage("gh:zfergus/filib#5b9cbc73790585b391da1624c1f1ea52eaaeb5ac")