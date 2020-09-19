# Copy header files into the build directory
function(ipc_toolkit_copy_headers)
  foreach(filepath IN ITEMS ${ARGN})
    file(RELATIVE_PATH filename "${CMAKE_SOURCE_DIR}/src" "${filepath}")
    if(${filename} MATCHES ".*\.(hpp|h|ipp|tpp)$")
      configure_file(${filepath} ${PROJECT_BINARY_DIR}/include/ipc/${filename})
    endif()
  endforeach()
endfunction()
