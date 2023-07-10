# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.27

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/andrew/Documents/GitHub/surface/HBV-BioRT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andrew/Documents/GitHub/surface/HBV-BioRT/build

# Include any dependencies generated for this target.
include dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/progress.make

# Include the compile flags for this target's objects.
include dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/flags.make

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/flags.make
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o: /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial.c
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o -MF CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o.d -o CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o -c /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial.c

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.i"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial.c > CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.i

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.s"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial.c -o CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.s

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/flags.make
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o: /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/test_sunlinsol.c
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o -MF CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o.d -o CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o -c /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/test_sunlinsol.c

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.i"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/test_sunlinsol.c > CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.i

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.s"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/test_sunlinsol.c -o CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.s

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/flags.make
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o: /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/src/sundials/sundials_linearsolver.c
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o -MF CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o.d -o CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o -c /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/src/sundials/sundials_linearsolver.c

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.i"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/src/sundials/sundials_linearsolver.c > CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.i

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.s"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/src/sundials/sundials_linearsolver.c -o CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.s

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/flags.make
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o: /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/src/sundials/sundials_nvector.c
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o -MF CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o.d -o CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o -c /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/src/sundials/sundials_nvector.c

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.i"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/src/sundials/sundials_nvector.c > CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.i

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.s"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/src/sundials/sundials_nvector.c -o CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.s

# Object files for target test_sunlinsol_sptfqmr_serial
test_sunlinsol_sptfqmr_serial_OBJECTS = \
"CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o" \
"CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o" \
"CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o" \
"CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o"

# External object files for target test_sunlinsol_sptfqmr_serial
test_sunlinsol_sptfqmr_serial_EXTERNAL_OBJECTS =

dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/test_sunlinsol_sptfqmr_serial.c.o
dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/test_sunlinsol.c.o
dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_linearsolver.c.o
dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/__/__/__/__/src/sundials/sundials_nvector.c.o
dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/build.make
dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: dep/cvode/src/nvector/serial/libsundials_nvecserial.so.4.1.0
dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: dep/cvode/src/sunlinsol/sptfqmr/libsundials_sunlinsolsptfqmr.so.2.1.0
dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: /usr/lib/x86_64-linux-gnu/librt.a
dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial: dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C executable test_sunlinsol_sptfqmr_serial"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/build: dep/cvode/examples/sunlinsol/sptfqmr/serial/test_sunlinsol_sptfqmr_serial
.PHONY : dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/build

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/clean:
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial && $(CMAKE_COMMAND) -P CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/cmake_clean.cmake
.PHONY : dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/clean

dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/depend:
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrew/Documents/GitHub/surface/HBV-BioRT /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/sunlinsol/sptfqmr/serial /home/andrew/Documents/GitHub/surface/HBV-BioRT/build /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : dep/cvode/examples/sunlinsol/sptfqmr/serial/CMakeFiles/test_sunlinsol_sptfqmr_serial.dir/depend

