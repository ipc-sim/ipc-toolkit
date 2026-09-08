# spdlog (https://github.com/gabime/spdlog)
# License: MIT
if(TARGET spdlog::spdlog)
    # Someone else created the target, so the fmt patch below never runs. That
    # is fine without CUDA, but the patch is what makes the bundled fmt
    # compile under nvcc at all, so warn rather than fail at the first .cu.
    if(IPC_TOOLKIT_WITH_CUDA)
        message(WARNING
            "spdlog::spdlog was provided by an enclosing project, so IPC "
            "Toolkit's cmake/patches/fmt-nvcc-compat.patch was not applied. "
            "The bundled fmt does not compile under nvcc unpatched: its "
            "literal-encoding probe misfires and a char32_t table uses hex "
            "escapes with the high bit set. Apply the same patch to your "
            "spdlog, or let IPC Toolkit fetch its own.")
    endif()
    return()
endif()

message(STATUS "Third-party: creating target 'spdlog::spdlog'")

option(SPDLOG_INSTALL "Generate the install target" ON)
set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME "spdlog")

include(CPM)
# The bundled fmt trips up NVCC's device front-end (EDG): a compile-time /utf-8
# probe misfires and a char32_t table uses hex escapes with the high bit set.
# Neither is fixable via compiler flags (the front-end never sees /utf-8), so we
# patch the bundled headers to guard both cases behind __CUDACC__.
CPMAddPackage(
    URI "gh:gabime/spdlog@1.17.0"
    PATCHES "${CMAKE_CURRENT_LIST_DIR}/../patches/fmt-nvcc-compat.patch"
)

set_target_properties(spdlog PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Folder name for IDE
set_target_properties(spdlog PROPERTIES FOLDER "ThirdParty")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" OR
   "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    target_compile_options(spdlog PRIVATE
        "-Wno-sign-conversion"
    )
endif()
