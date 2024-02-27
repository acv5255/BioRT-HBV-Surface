cmake_minimum_required(VERSION 3.22)

project(
    BioRT
    LANGUAGES C
)

if (WIN32)
    add_compile_options(/MT /w)
endif()

add_subdirectory(dep/cvode)
add_subdirectory(src)
