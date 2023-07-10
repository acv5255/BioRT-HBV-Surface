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
include dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/progress.make

# Include the compile flags for this target's objects.
include dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/flags.make

dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o: dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/flags.make
dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o: /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/cvode/serial/cvDisc_dns.c
dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o: dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/cvode/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o -MF CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o.d -o CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o -c /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/cvode/serial/cvDisc_dns.c

dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.i"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/cvode/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/cvode/serial/cvDisc_dns.c > CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.i

dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.s"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/cvode/serial && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/cvode/serial/cvDisc_dns.c -o CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.s

# Object files for target cvDisc_dns
cvDisc_dns_OBJECTS = \
"CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o"

# External object files for target cvDisc_dns
cvDisc_dns_EXTERNAL_OBJECTS =

dep/cvode/examples/cvode/serial/cvDisc_dns: dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/cvDisc_dns.c.o
dep/cvode/examples/cvode/serial/cvDisc_dns: dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/build.make
dep/cvode/examples/cvode/serial/cvDisc_dns: dep/cvode/src/cvode/libsundials_cvode.so.4.1.0
dep/cvode/examples/cvode/serial/cvDisc_dns: dep/cvode/src/nvector/serial/libsundials_nvecserial.so.4.1.0
dep/cvode/examples/cvode/serial/cvDisc_dns: /usr/lib/x86_64-linux-gnu/librt.a
dep/cvode/examples/cvode/serial/cvDisc_dns: dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/andrew/Documents/GitHub/surface/HBV-BioRT/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable cvDisc_dns"
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/cvode/serial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cvDisc_dns.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/build: dep/cvode/examples/cvode/serial/cvDisc_dns
.PHONY : dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/build

dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/clean:
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/cvode/serial && $(CMAKE_COMMAND) -P CMakeFiles/cvDisc_dns.dir/cmake_clean.cmake
.PHONY : dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/clean

dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/depend:
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrew/Documents/GitHub/surface/HBV-BioRT /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode/examples/cvode/serial /home/andrew/Documents/GitHub/surface/HBV-BioRT/build /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/cvode/serial /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : dep/cvode/examples/cvode/serial/CMakeFiles/cvDisc_dns.dir/depend

