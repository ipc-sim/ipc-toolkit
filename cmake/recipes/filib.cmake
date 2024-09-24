# filib (https://github.com/zfergus/filib.git)
# License: LGPL-2.1
if(TARGET filib::filib)
  return()
endif()

message(STATUS "Third-party: creating target 'filib::filib'")

# filib should be built as a shared library to avoid licensing IPC Toolkit under LGPL
if(WIN32 AND NOT IPC_TOOLKIT_TOPLEVEL_PROJECT)
  # Setting up proper linkage on Windows is a bit tricky, so we'll just use a
  # static library by default and provide instructions on how to build as a
  # shared library in the README.
  option(FILIB_BUILD_SHARED_LIB "Build shared library" OFF)
else()
  # NOTE: Our Windows CMake is properly configured to build shared libraries for
  #       OUR applications and python bindings.
  option(FILIB_BUILD_SHARED_LIB "Build shared library" ON)
endif()

include(CPM)
CPMAddPackage("gh:zfergus/filib#03e4eb0fc59399bd0003f8efd3179078195df49f")