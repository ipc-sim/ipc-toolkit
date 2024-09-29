# filib (https://github.com/zfergus/filib.git)
# License: LGPL-2.1
if(TARGET filib::filib)
  return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

# NOTE: filib should be built as a shared library to avoid licensing IPC Toolkit under LGPL
# However, Setting up proper linkage is a bit tricky, so we'll just use a static library by
# default and provide instructions on how to build as a shared library in the README.
option(FILIB_BUILD_SHARED_LIB "Build shared library" OFF)

include(CPM)
CPMAddPackage("gh:zfergus/filib#03e4eb0fc59399bd0003f8efd3179078195df49f")