###############################################################################
# Create Python module for IPC Toolkit
include(pybind11)

pybind11_add_module(ipctk)

add_subdirectory(src) # Add sources to ipctk
target_include_directories(ipctk PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

target_link_libraries(ipctk PRIVATE ipc::toolkit)

# libigl (binding igl::edges)
include(libigl)
target_link_libraries(ipctk PRIVATE igl::core igl::predicates)

# Tight Inclusion CCD (binding tight_inclusion as ipctk.tight_inclusion)
include(tight_inclusion)
target_link_libraries(ipctk PRIVATE tight_inclusion::tight_inclusion)

# TBB (binding tbb::global_control)
include(onetbb)
target_link_libraries(ipctk PRIVATE TBB::tbb)

if (FILIB_BUILD_SHARED_LIB AND WIN32)
  # Copy DLLs to the output directory
  add_custom_command(
    TARGET ipctk POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_RUNTIME_DLLS:ipctk> $<TARGET_FILE_DIR:ipctk>
    COMMAND_EXPAND_LISTS
  )
endif()

if (IPC_TOOLKIT_WITH_PROFILER)
  include(pybind11_json)
  target_link_libraries(ipctk PRIVATE pybind11::json)
endif()

# Extra warnings
# target_link_libraries(ipctk PRIVATE IPCToolkit::warnings)

# Folder name for IDE
set_target_properties(ipctk PROPERTIES FOLDER "SRC")
get_target_property(IPCTK_SOURCES ipctk SOURCES)
source_group(TREE "${PROJECT_SOURCE_DIR}/python" FILES ${IPCTK_SOURCES})