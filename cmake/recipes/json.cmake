# Nlohmann's JSON (https://github.com/nlohmann/json)
# License: MIT
if(TARGET nlohmann_json::nlohmann_json)
    return()
endif()

message(STATUS "Third-party: creating target 'nlohmann_json::nlohmann_json'")

# nlohmann_json is a big repo for a single header, so we just download the release archive
set(NLOHMANNJSON_VERSION "v3.11.2")

include(CPM)
CPMAddPackage(
    NAME nlohmann_json
    URL "https://github.com/nlohmann/json/releases/download/${NLOHMANNJSON_VERSION}/include.zip"
    URL_HASH SHA256=e5c7a9f49a16814be27e4ed0ee900ecd0092bfb7dbfca65b5a421b774dccaaed
    DOWNLOAD_ONLY YES
)

add_library(nlohmann_json INTERFACE)
add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)

include(GNUInstallDirs)
target_include_directories(nlohmann_json SYSTEM INTERFACE
    "$<BUILD_INTERFACE:${nlohmann_json_SOURCE_DIR}>/include"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Install rules
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME nlohmann_json)
install(DIRECTORY ${nlohmann_json_SOURCE_DIR}/include/nlohmann DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS nlohmann_json EXPORT NlohmannJson_Targets)
install(EXPORT NlohmannJson_Targets DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nlohmann_json NAMESPACE nlohmann_json::)
export(EXPORT NlohmannJson_Targets FILE "${CMAKE_CURRENT_BINARY_DIR}/NlohmannJsonTargets.cmake")