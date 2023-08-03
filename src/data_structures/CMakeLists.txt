cmake_minimum_required(VERSION 3.26)

add_library(data_structures
    CalibrationStruct.cpp    
    ChemicalState.cpp
    ChemTableEntry.cpp
    ControlData.cpp
    HBVParameters.cpp
    KineticTableEntry.cpp
    ReactionNetwork.cpp
    SoilConstants.cpp
    SoilParameters.cpp
    Subcatchment.cpp
)

set_property(TARGET data_structures PROPERTY CXX_STANDARD 20)

target_compile_options(
    data_structures
    PRIVATE
    -Wall
    -Wextra
    -Werror
)

target_include_directories(
    data_structures PRIVATE
    ${BioRT_SOURCE_DIR}/build/dep/cvode/include
    ${BioRT_SOURCE_DIR}/dep/cvode/include
    ${BioRT_SOURCE_DIR}/include
)
