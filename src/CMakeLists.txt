cmake_minimum_required(VERSION 3.26)

project(biort)

add_subdirectory(data_structures)

add_executable(biort
    main.cpp
    custom_io.cpp
    init.cpp
    lookup.cpp
    optparse.cpp
    print.cpp
    react.cpp
    read_chem.cpp
    read_cini.cpp
    read_hbv.cpp
    read_param.cpp
    read_precipchem.cpp
    read_soil.cpp
    set_numexp.cpp
    speciation.cpp
    time_func.cpp
    transport.cpp
    util_func.cpp
)

set_property(TARGET biort PROPERTY CXX_STANDARD 20)

target_compile_options(
    biort
    PRIVATE
    -Wall
    -Wextra
    -Werror
)

target_include_directories(
    biort PRIVATE
    ${BioRT_SOURCE_DIR}/build/dep/cvode/include
    ${BioRT_SOURCE_DIR}/dep/cvode/include
    ${BioRT_SOURCE_DIR}/include
)

target_link_directories(
    biort PRIVATE
    ${BioRT_SOURCE_DIR}/build/dep/cvode/src/cvode
    ${BioRT_SOURCE_DIR}/build/dep/cvode/src/nvector/serial
)

target_link_libraries(
    biort PRIVATE
    m
    sundials_cvode
    sundials_nvecserial
    data_structures
)