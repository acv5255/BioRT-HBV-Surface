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
include dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/progress.make

# Include the compile flags for this target's objects.
include dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/flags.make

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/flags.make
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse.c
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o -MF CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o.d -o CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse.c

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse.c > CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.i

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse.c -o CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.s

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/flags.make
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/test_sunmatrix.c
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o -MF CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o.d -o CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/test_sunmatrix.c

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/test_sunmatrix.c > CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.i

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/test_sunmatrix.c -o CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.s

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/flags.make
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_matrix.c
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o -MF CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o.d -o CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_matrix.c

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_matrix.c > CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.i

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_matrix.c -o CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.s

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/flags.make
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_nvector.c
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o -MF CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o.d -o CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_nvector.c

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_nvector.c > CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.i

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_nvector.c -o CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.s

# Object files for target test_sunmatrix_sparse
test_sunmatrix_sparse_OBJECTS = \
"CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o" \
"CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o" \
"CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o" \
"CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o"

# External object files for target test_sunmatrix_sparse
test_sunmatrix_sparse_EXTERNAL_OBJECTS =

dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/test_sunmatrix_sparse.c.o
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/test_sunmatrix.c.o
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_matrix.c.o
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/__/__/__/src/sundials/sundials_nvector.c.o
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/build.make
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/src/nvector/serial/libsundials_nvecserial.so.4.1.0
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/src/sunmatrix/dense/libsundials_sunmatrixdense.so.2.1.0
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/src/sunmatrix/band/libsundials_sunmatrixband.so.2.1.0
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/src/sunmatrix/sparse/libsundials_sunmatrixsparse.so.2.1.0
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: /usr/lib/x86_64-linux-gnu/librt.a
dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse: dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C executable test_sunmatrix_sparse"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test_sunmatrix_sparse.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/build: dep/cvode/examples/sunmatrix/sparse/test_sunmatrix_sparse
.PHONY : dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/build

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/clean:
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse && $(CMAKE_COMMAND) -P CMakeFiles/test_sunmatrix_sparse.dir/cmake_clean.cmake
.PHONY : dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/clean

dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/depend:
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/examples/sunmatrix/sparse /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : dep/cvode/examples/sunmatrix/sparse/CMakeFiles/test_sunmatrix_sparse.dir/depend

