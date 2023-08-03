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
include dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/progress.make

# Include the compile flags for this target's objects.
include dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/flags.make

dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o: dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/flags.make
dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol.c
dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o: dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/cvode/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o -MF CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o.d -o CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol.c

dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/cvode/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol.c > CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.i

dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/cvode/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol.c -o CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.s

# Object files for target cvRoberts_dns_negsol
cvRoberts_dns_negsol_OBJECTS = \
"CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o"

# External object files for target cvRoberts_dns_negsol
cvRoberts_dns_negsol_EXTERNAL_OBJECTS =

dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol: dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/cvRoberts_dns_negsol.c.o
dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol: dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/build.make
dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol: dep/cvode/src/cvode/libsundials_cvode.so.4.1.0
dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol: dep/cvode/src/nvector/serial/libsundials_nvecserial.so.4.1.0
dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol: /usr/lib/x86_64-linux-gnu/librt.a
dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol: dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable cvRoberts_dns_negsol"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/cvode/serial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cvRoberts_dns_negsol.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/build: dep/cvode/examples/cvode/serial/cvRoberts_dns_negsol
.PHONY : dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/build

dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/clean:
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/cvode/serial && $(CMAKE_COMMAND) -P CMakeFiles/cvRoberts_dns_negsol.dir/cmake_clean.cmake
.PHONY : dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/clean

dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/depend:
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/cvode/serial /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/cvode/serial /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : dep/cvode/examples/cvode/serial/CMakeFiles/cvRoberts_dns_negsol.dir/depend
