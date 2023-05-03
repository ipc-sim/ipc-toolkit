if(TARGET broadmark)
    return()
endif()

message(STATUS "Third-party: creating target 'broadmark'")

include(FetchContent)
FetchContent_Declare(
    broadmark
    GIT_REPOSITORY https://github.com/dbelgrod/Broadmark.git
    GIT_TAG master
    GIT_SHALLOW FALSE
)
FetchContent_GetProperties(broadmark)

if(NOT broadmark_POPULATED)
    FetchContent_Populate(broadmark)
    add_subdirectory(${broadmark_SOURCE_DIR} ${broadmark_BINARY_DIR})
endif()

file(GLOB BROADMARK_SOURCES
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/*.cpp"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/*.cpp"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/*.cpp"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/Bullet3Common/*.cpp"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/Bullet3Collision/BroadPhaseCollision/*.cpp"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet2/*.cpp"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet2/LinearMath/*.cpp"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet2/BulletCollision/*.cpp"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet2/BulletDynamics/*.cpp"
)

# Create an imported target for the Broadmark library
add_library(broadmark STATIC ${BROADMARK_SOURCES})
target_include_directories(broadmark PUBLIC
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies" # Add Bullet3 include directory
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/Bullet3Common"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet3/Bullet3Collision/BroadPhaseCollision"
    "${broadmark_SOURCE_DIR}/Broadmark/Algorithms/Dependencies/Bullet2"
)
