# pybind11_json (https://github.com/pybind/pybind11_json)
# License: MIT
if(TARGET pybind11::json)
    return()
endif()

message(STATUS "Third-party: creating target 'pybind11::json'")

include(CPM)
CPMAddPackage(
    URI "gh:pybind/pybind11_json#0.2.15"
    DOWNLOAD_ONLY YES
)

add_library(pybind11_json INTERFACE)
add_library(pybind11::json ALIAS pybind11_json)

target_include_directories(pybind11_json INTERFACE
    "$<BUILD_INTERFACE:${pybind11_json_SOURCE_DIR}/include>"
)

include(pybind11)
target_link_libraries(pybind11_json INTERFACE pybind11::pybind11)

include(json)
target_link_libraries(pybind11_json INTERFACE nlohmann_json::nlohmann_json)