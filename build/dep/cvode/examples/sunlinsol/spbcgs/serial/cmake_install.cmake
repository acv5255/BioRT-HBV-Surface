# Install script for directory: /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/spbcgs/serial/test_sunlinsol_spbcgs_serial.c;/usr/local/examples/sunlinsol/spbcgs/serial/test_sunlinsol.h;/usr/local/examples/sunlinsol/spbcgs/serial/test_sunlinsol.c;/usr/local/examples/sunlinsol/spbcgs/serial/sundials_linearsolver.c;/usr/local/examples/sunlinsol/spbcgs/serial/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/spbcgs/serial" TYPE FILE FILES
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/test_sunlinsol_spbcgs_serial.c"
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/../../test_sunlinsol.h"
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/../../test_sunlinsol.c"
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/../../../../src/sundials/sundials_linearsolver.c"
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/../../../../src/sundials/sundials_nvector.c"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/spbcgs/serial/test_sunlinsol_spbcgs_serial.c;/usr/local/examples/sunlinsol/spbcgs/serial/test_sunlinsol.h;/usr/local/examples/sunlinsol/spbcgs/serial/test_sunlinsol.c;/usr/local/examples/sunlinsol/spbcgs/serial/sundials_linearsolver.c;/usr/local/examples/sunlinsol/spbcgs/serial/sundials_nvector.c")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/spbcgs/serial" TYPE FILE FILES
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/test_sunlinsol_spbcgs_serial.c"
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/../../test_sunlinsol.h"
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/../../test_sunlinsol.c"
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/../../../../src/sundials/sundials_linearsolver.c"
    "/home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/spbcgs/serial/../../../../src/sundials/sundials_nvector.c"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/spbcgs/serial/CMakeLists.txt")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/spbcgs/serial" TYPE FILE FILES "/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/spbcgs/serial/CMakeLists.txt")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/examples/sunlinsol/spbcgs/serial/Makefile")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/examples/sunlinsol/spbcgs/serial" TYPE FILE RENAME "Makefile" FILES "/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/spbcgs/serial/Makefile_ex")
endif()

