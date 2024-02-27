cmake_minimum_required(VERSION 3.22)

project(biort)

add_executable(biort
    main.c
    custom_io.c
    init.c
    lookup.c
    optparse.c
    print.c
    react.c
    read_chem.c
    read_cini.c
    read_hbv.c
    read_param.c
    read_precipchem.c
    read_soil.c
    set_numexp.c
    speciation.c
    time_func.c
    transpt.c
    util_func.c
)


IF(WIN32)
    add_compile_definitions(
        /DWIN32     
    )
ELSE()
    add_compile_definitions(
        O2
    )
ENDIF()

if(WIN32)
    add_compile_options(
        /MT
    )
endif()


target_compile_options(
    biort
    PRIVATE
    -Wall
    -Wextra
    # -Werror
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

# target_link_libraries(
#     biort PRIVATE
#     m
#     sundials_cvode
#     sundials_nvecserial
# )


if (WIN32)
	target_link_libraries(
	    biort PRIVATE
        # libucrt
	    libvcruntime
	    sundials_cvode
	    sundials_nvecserial
	)
else()
	target_link_libraries(
	    biort PRIVATE
        m
	    sundials_cvode
	    sundials_nvecserial
	)
endif()

add_dependencies(biort sundials_cvode_static)
add_dependencies(biort sundials_nvecserial_static)

add_custom_command(
    TARGET biort
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:biort> ${BioRT_SOURCE_DIR}
)