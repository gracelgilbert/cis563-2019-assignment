# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.15.3/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.15.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build

# Include any dependencies generated for this target.
include Projects/mass_spring/CMakeFiles/mass_spring.dir/depend.make

# Include the progress variables for this target.
include Projects/mass_spring/CMakeFiles/mass_spring.dir/progress.make

# Include the compile flags for this target's objects.
include Projects/mass_spring/CMakeFiles/mass_spring.dir/flags.make

Projects/mass_spring/CMakeFiles/mass_spring.dir/main.cpp.o: Projects/mass_spring/CMakeFiles/mass_spring.dir/flags.make
Projects/mass_spring/CMakeFiles/mass_spring.dir/main.cpp.o: ../Projects/mass_spring/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Projects/mass_spring/CMakeFiles/mass_spring.dir/main.cpp.o"
	cd /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/Projects/mass_spring && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/mass_spring.dir/main.cpp.o -c /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/Projects/mass_spring/main.cpp

Projects/mass_spring/CMakeFiles/mass_spring.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mass_spring.dir/main.cpp.i"
	cd /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/Projects/mass_spring && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/Projects/mass_spring/main.cpp > CMakeFiles/mass_spring.dir/main.cpp.i

Projects/mass_spring/CMakeFiles/mass_spring.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mass_spring.dir/main.cpp.s"
	cd /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/Projects/mass_spring && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/Projects/mass_spring/main.cpp -o CMakeFiles/mass_spring.dir/main.cpp.s

# Object files for target mass_spring
mass_spring_OBJECTS = \
"CMakeFiles/mass_spring.dir/main.cpp.o"

# External object files for target mass_spring
mass_spring_EXTERNAL_OBJECTS =

../Projects/mass_spring/mass_spring: Projects/mass_spring/CMakeFiles/mass_spring.dir/main.cpp.o
../Projects/mass_spring/mass_spring: Projects/mass_spring/CMakeFiles/mass_spring.dir/build.make
../Projects/mass_spring/mass_spring: Projects/mass_spring/CMakeFiles/mass_spring.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../../Projects/mass_spring/mass_spring"
	cd /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/Projects/mass_spring && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mass_spring.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Projects/mass_spring/CMakeFiles/mass_spring.dir/build: ../Projects/mass_spring/mass_spring

.PHONY : Projects/mass_spring/CMakeFiles/mass_spring.dir/build

Projects/mass_spring/CMakeFiles/mass_spring.dir/clean:
	cd /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/Projects/mass_spring && $(CMAKE_COMMAND) -P CMakeFiles/mass_spring.dir/cmake_clean.cmake
.PHONY : Projects/mass_spring/CMakeFiles/mass_spring.dir/clean

Projects/mass_spring/CMakeFiles/mass_spring.dir/depend:
	cd /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/Projects/mass_spring /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/Projects/mass_spring /Users/gracegilbert/Documents/CIS563/cis563-2019-assignment/build/Projects/mass_spring/CMakeFiles/mass_spring.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Projects/mass_spring/CMakeFiles/mass_spring.dir/depend

