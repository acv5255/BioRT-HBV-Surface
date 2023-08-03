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
CMAKE_SOURCE_DIR = /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build

# Include any dependencies generated for this target.
include dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/progress.make

# Include the compile flags for this target's objects.
include dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/flags.make

dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o: dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/flags.make
dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/nvector/serial/nvector_serial.c
dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o: dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o -MF CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o.d -o CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/nvector/serial/nvector_serial.c

dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/nvector/serial/nvector_serial.c > CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.i

dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/nvector/serial/nvector_serial.c -o CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.s

dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o: dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/flags.make
dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_math.c
dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o: dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o -MF CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o.d -o CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_math.c

dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_math.c > CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.i

dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_math.c -o CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.s

# Object files for target sundials_nvecserial_static
sundials_nvecserial_static_OBJECTS = \
"CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o" \
"CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o"

# External object files for target sundials_nvecserial_static
sundials_nvecserial_static_EXTERNAL_OBJECTS =

dep/cvode/src/nvector/serial/libsundials_nvecserial.a: dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/nvector_serial.c.o
dep/cvode/src/nvector/serial/libsundials_nvecserial.a: dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/__/__/sundials/sundials_math.c.o
dep/cvode/src/nvector/serial/libsundials_nvecserial.a: dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/build.make
dep/cvode/src/nvector/serial/libsundials_nvecserial.a: dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C static library libsundials_nvecserial.a"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && $(CMAKE_COMMAND) -P CMakeFiles/sundials_nvecserial_static.dir/cmake_clean_target.cmake
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sundials_nvecserial_static.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/build: dep/cvode/src/nvector/serial/libsundials_nvecserial.a
.PHONY : dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/build

dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/clean:
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial && $(CMAKE_COMMAND) -P CMakeFiles/sundials_nvecserial_static.dir/cmake_clean.cmake
.PHONY : dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/clean

dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/depend:
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/nvector/serial /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : dep/cvode/src/nvector/serial/CMakeFiles/sundials_nvecserial_static.dir/depend

