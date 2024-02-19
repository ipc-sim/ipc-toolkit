###############################################################################
# Create Python module for IPC Toolkit
include(pybind11)

pybind11_add_module(ipctk)

add_subdirectory(src) # Add sources to ipctk
target_include_directories(ipctk PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(ipctk PRIVATE ipc::toolkit)

# Extra warnings
# target_link_libraries(ipctk PRIVATE IPCToolkit::warnings)