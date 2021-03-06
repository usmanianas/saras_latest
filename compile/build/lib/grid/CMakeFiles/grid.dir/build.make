# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /home/anas/local/bin/cmake

# The command to remove a file.
RM = /home/anas/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build

# Include any dependencies generated for this target.
include lib/grid/CMakeFiles/grid.dir/depend.make

# Include the progress variables for this target.
include lib/grid/CMakeFiles/grid.dir/progress.make

# Include the compile flags for this target's objects.
include lib/grid/CMakeFiles/grid.dir/flags.make

lib/grid/CMakeFiles/grid.dir/grid.cc.o: lib/grid/CMakeFiles/grid.dir/flags.make
lib/grid/CMakeFiles/grid.dir/grid.cc.o: ../../lib/grid/grid.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/grid/CMakeFiles/grid.dir/grid.cc.o"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/grid && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/grid.dir/grid.cc.o -c /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/grid/grid.cc

lib/grid/CMakeFiles/grid.dir/grid.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/grid.dir/grid.cc.i"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/grid && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/grid/grid.cc > CMakeFiles/grid.dir/grid.cc.i

lib/grid/CMakeFiles/grid.dir/grid.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/grid.dir/grid.cc.s"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/grid && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/grid/grid.cc -o CMakeFiles/grid.dir/grid.cc.s

lib/grid/CMakeFiles/grid.dir/grid.cc.o.requires:

.PHONY : lib/grid/CMakeFiles/grid.dir/grid.cc.o.requires

lib/grid/CMakeFiles/grid.dir/grid.cc.o.provides: lib/grid/CMakeFiles/grid.dir/grid.cc.o.requires
	$(MAKE) -f lib/grid/CMakeFiles/grid.dir/build.make lib/grid/CMakeFiles/grid.dir/grid.cc.o.provides.build
.PHONY : lib/grid/CMakeFiles/grid.dir/grid.cc.o.provides

lib/grid/CMakeFiles/grid.dir/grid.cc.o.provides.build: lib/grid/CMakeFiles/grid.dir/grid.cc.o


# Object files for target grid
grid_OBJECTS = \
"CMakeFiles/grid.dir/grid.cc.o"

# External object files for target grid
grid_EXTERNAL_OBJECTS =

lib/grid/libgrid.a: lib/grid/CMakeFiles/grid.dir/grid.cc.o
lib/grid/libgrid.a: lib/grid/CMakeFiles/grid.dir/build.make
lib/grid/libgrid.a: lib/grid/CMakeFiles/grid.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libgrid.a"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/grid && $(CMAKE_COMMAND) -P CMakeFiles/grid.dir/cmake_clean_target.cmake
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/grid && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/grid.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/grid/CMakeFiles/grid.dir/build: lib/grid/libgrid.a

.PHONY : lib/grid/CMakeFiles/grid.dir/build

lib/grid/CMakeFiles/grid.dir/requires: lib/grid/CMakeFiles/grid.dir/grid.cc.o.requires

.PHONY : lib/grid/CMakeFiles/grid.dir/requires

lib/grid/CMakeFiles/grid.dir/clean:
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/grid && $(CMAKE_COMMAND) -P CMakeFiles/grid.dir/cmake_clean.cmake
.PHONY : lib/grid/CMakeFiles/grid.dir/clean

lib/grid/CMakeFiles/grid.dir/depend:
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/grid /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/grid /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/grid/CMakeFiles/grid.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/grid/CMakeFiles/grid.dir/depend

