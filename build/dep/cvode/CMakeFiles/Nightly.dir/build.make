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

# Utility rule file for Nightly.

# Include any custom commands dependencies for this target.
include dep/cvode/CMakeFiles/Nightly.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/cvode/CMakeFiles/Nightly.dir/progress.make

dep/cvode/CMakeFiles/Nightly:
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode && /usr/local/bin/ctest -D Nightly

Nightly: dep/cvode/CMakeFiles/Nightly
Nightly: dep/cvode/CMakeFiles/Nightly.dir/build.make
.PHONY : Nightly

# Rule to build all files generated by this target.
dep/cvode/CMakeFiles/Nightly.dir/build: Nightly
.PHONY : dep/cvode/CMakeFiles/Nightly.dir/build

dep/cvode/CMakeFiles/Nightly.dir/clean:
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode && $(CMAKE_COMMAND) -P CMakeFiles/Nightly.dir/cmake_clean.cmake
.PHONY : dep/cvode/CMakeFiles/Nightly.dir/clean

dep/cvode/CMakeFiles/Nightly.dir/depend:
	cd /home/andrew/Documents/GitHub/surface/HBV-BioRT/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrew/Documents/GitHub/surface/HBV-BioRT /home/andrew/Documents/GitHub/surface/HBV-BioRT/dep/cvode /home/andrew/Documents/GitHub/surface/HBV-BioRT/build /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode /home/andrew/Documents/GitHub/surface/HBV-BioRT/build/dep/cvode/CMakeFiles/Nightly.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : dep/cvode/CMakeFiles/Nightly.dir/depend

