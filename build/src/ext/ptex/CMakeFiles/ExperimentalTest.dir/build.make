# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build

# Utility rule file for ExperimentalTest.

# Include the progress variables for this target.
include src/ext/ptex/CMakeFiles/ExperimentalTest.dir/progress.make

src/ext/ptex/CMakeFiles/ExperimentalTest:
	cd /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex && /usr/bin/ctest -D ExperimentalTest

ExperimentalTest: src/ext/ptex/CMakeFiles/ExperimentalTest
ExperimentalTest: src/ext/ptex/CMakeFiles/ExperimentalTest.dir/build.make

.PHONY : ExperimentalTest

# Rule to build all files generated by this target.
src/ext/ptex/CMakeFiles/ExperimentalTest.dir/build: ExperimentalTest

.PHONY : src/ext/ptex/CMakeFiles/ExperimentalTest.dir/build

src/ext/ptex/CMakeFiles/ExperimentalTest.dir/clean:
	cd /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex && $(CMAKE_COMMAND) -P CMakeFiles/ExperimentalTest.dir/cmake_clean.cmake
.PHONY : src/ext/ptex/CMakeFiles/ExperimentalTest.dir/clean

src/ext/ptex/CMakeFiles/ExperimentalTest.dir/depend:
	cd /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3 /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/src/ext/ptex /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex /home/dennis/Documents/phd/courses/topics-in-computer-graphics/project/code/repo/pbrt-v3/build/src/ext/ptex/CMakeFiles/ExperimentalTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/ext/ptex/CMakeFiles/ExperimentalTest.dir/depend

