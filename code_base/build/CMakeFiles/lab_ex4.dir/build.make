# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_COMMAND = /usr/local/lib/python3.8/dist-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /usr/local/lib/python3.8/dist-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/e/GitHub/UdeS_S6_APP1/code_base

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/e/GitHub/UdeS_S6_APP1/code_base/build

# Include any dependencies generated for this target.
include CMakeFiles/lab_ex4.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/lab_ex4.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/lab_ex4.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lab_ex4.dir/flags.make

CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o: CMakeFiles/lab_ex4.dir/flags.make
CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o: /mnt/e/GitHub/UdeS_S6_APP1/code_base/src/lab_ex4.cpp
CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o: CMakeFiles/lab_ex4.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/e/GitHub/UdeS_S6_APP1/code_base/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o"
	/usr/lib/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o -MF CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o.d -o CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o -c /mnt/e/GitHub/UdeS_S6_APP1/code_base/src/lab_ex4.cpp

CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.i"
	/usr/lib/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/e/GitHub/UdeS_S6_APP1/code_base/src/lab_ex4.cpp > CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.i

CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.s"
	/usr/lib/ccache/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/e/GitHub/UdeS_S6_APP1/code_base/src/lab_ex4.cpp -o CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.s

# Object files for target lab_ex4
lab_ex4_OBJECTS = \
"CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o"

# External object files for target lab_ex4
lab_ex4_EXTERNAL_OBJECTS =

lab_ex4: CMakeFiles/lab_ex4.dir/src/lab_ex4.cpp.o
lab_ex4: CMakeFiles/lab_ex4.dir/build.make
lab_ex4: CMakeFiles/lab_ex4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/e/GitHub/UdeS_S6_APP1/code_base/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable lab_ex4"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lab_ex4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lab_ex4.dir/build: lab_ex4
.PHONY : CMakeFiles/lab_ex4.dir/build

CMakeFiles/lab_ex4.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lab_ex4.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lab_ex4.dir/clean

CMakeFiles/lab_ex4.dir/depend:
	cd /mnt/e/GitHub/UdeS_S6_APP1/code_base/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/e/GitHub/UdeS_S6_APP1/code_base /mnt/e/GitHub/UdeS_S6_APP1/code_base /mnt/e/GitHub/UdeS_S6_APP1/code_base/build /mnt/e/GitHub/UdeS_S6_APP1/code_base/build /mnt/e/GitHub/UdeS_S6_APP1/code_base/build/CMakeFiles/lab_ex4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lab_ex4.dir/depend

