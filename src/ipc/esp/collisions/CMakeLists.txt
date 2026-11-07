set(SOURCES
    esp_collision.cpp
    esp_collision.hpp
    vertex_matrix_view.hpp
    esp_collision_template.cpp
    esp_collision_template.hpp
    esp_collision_dict.cpp
    esp_collision_dict.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ipc_toolkit PRIVATE ${SOURCES})
