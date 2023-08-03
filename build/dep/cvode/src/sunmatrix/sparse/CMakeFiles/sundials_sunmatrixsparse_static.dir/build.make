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
include dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/compiler_depend.make

# Include the progress variables for this target.
include dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/progress.make

# Include the compile flags for this target's objects.
include dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/flags.make

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/flags.make
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sunmatrix/sparse/sunmatrix_sparse.c
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o -MF CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o.d -o CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sunmatrix/sparse/sunmatrix_sparse.c

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sunmatrix/sparse/sunmatrix_sparse.c > CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.i

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sunmatrix/sparse/sunmatrix_sparse.c -o CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.s

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/flags.make
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_nvector.c
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o -MF CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o.d -o CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_nvector.c

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_nvector.c > CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.i

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_nvector.c -o CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.s

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/flags.make
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_matrix.c
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o -MF CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o.d -o CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_matrix.c

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_matrix.c > CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.i

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_matrix.c -o CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.s

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/flags.make
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o: /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_math.c
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o -MF CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o.d -o CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o -c /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_math.c

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing C source to CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.i"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_math.c > CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.i

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling C source to assembly CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.s"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && /usr/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sundials/sundials_math.c -o CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.s

# Object files for target sundials_sunmatrixsparse_static
sundials_sunmatrixsparse_static_OBJECTS = \
"CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o" \
"CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o" \
"CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o" \
"CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o"

# External object files for target sundials_sunmatrixsparse_static
sundials_sunmatrixsparse_static_EXTERNAL_OBJECTS =

dep/cvode/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/sunmatrix_sparse.c.o
dep/cvode/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_nvector.c.o
dep/cvode/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_matrix.c.o
dep/cvode/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/__/__/sundials/sundials_math.c.o
dep/cvode/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/build.make
dep/cvode/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a: dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking C static library libsundials_sunmatrixsparse.a"
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && $(CMAKE_COMMAND) -P CMakeFiles/sundials_sunmatrixsparse_static.dir/cmake_clean_target.cmake
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sundials_sunmatrixsparse_static.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/build: dep/cvode/src/sunmatrix/sparse/libsundials_sunmatrixsparse.a
.PHONY : dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/build

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/clean:
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse && $(CMAKE_COMMAND) -P CMakeFiles/sundials_sunmatrixsparse_static.dir/cmake_clean.cmake
.PHONY : dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/clean

dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/depend:
	cd /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/dep/cvode/src/sunmatrix/sparse /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse /home/andrew/Documents/GitHub/biort-surface-cpp/BioRT-HBV-Surface/build/dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : dep/cvode/src/sunmatrix/sparse/CMakeFiles/sundials_sunmatrixsparse_static.dir/depend
