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
include lib/io/CMakeFiles/tseries.dir/depend.make

# Include the progress variables for this target.
include lib/io/CMakeFiles/tseries.dir/progress.make

# Include the compile flags for this target's objects.
include lib/io/CMakeFiles/tseries.dir/flags.make

lib/io/CMakeFiles/tseries.dir/tseries.cc.o: lib/io/CMakeFiles/tseries.dir/flags.make
lib/io/CMakeFiles/tseries.dir/tseries.cc.o: ../../lib/io/tseries.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/io/CMakeFiles/tseries.dir/tseries.cc.o"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/io && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tseries.dir/tseries.cc.o -c /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/io/tseries.cc

lib/io/CMakeFiles/tseries.dir/tseries.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tseries.dir/tseries.cc.i"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/io && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/io/tseries.cc > CMakeFiles/tseries.dir/tseries.cc.i

lib/io/CMakeFiles/tseries.dir/tseries.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tseries.dir/tseries.cc.s"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/io && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/io/tseries.cc -o CMakeFiles/tseries.dir/tseries.cc.s

lib/io/CMakeFiles/tseries.dir/tseries.cc.o.requires:

.PHONY : lib/io/CMakeFiles/tseries.dir/tseries.cc.o.requires

lib/io/CMakeFiles/tseries.dir/tseries.cc.o.provides: lib/io/CMakeFiles/tseries.dir/tseries.cc.o.requires
	$(MAKE) -f lib/io/CMakeFiles/tseries.dir/build.make lib/io/CMakeFiles/tseries.dir/tseries.cc.o.provides.build
.PHONY : lib/io/CMakeFiles/tseries.dir/tseries.cc.o.provides

lib/io/CMakeFiles/tseries.dir/tseries.cc.o.provides.build: lib/io/CMakeFiles/tseries.dir/tseries.cc.o


# Object files for target tseries
tseries_OBJECTS = \
"CMakeFiles/tseries.dir/tseries.cc.o"

# External object files for target tseries
tseries_EXTERNAL_OBJECTS =

lib/io/libtseries.a: lib/io/CMakeFiles/tseries.dir/tseries.cc.o
lib/io/libtseries.a: lib/io/CMakeFiles/tseries.dir/build.make
lib/io/libtseries.a: lib/io/CMakeFiles/tseries.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libtseries.a"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/io && $(CMAKE_COMMAND) -P CMakeFiles/tseries.dir/cmake_clean_target.cmake
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/io && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tseries.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/io/CMakeFiles/tseries.dir/build: lib/io/libtseries.a

.PHONY : lib/io/CMakeFiles/tseries.dir/build

lib/io/CMakeFiles/tseries.dir/requires: lib/io/CMakeFiles/tseries.dir/tseries.cc.o.requires

.PHONY : lib/io/CMakeFiles/tseries.dir/requires

lib/io/CMakeFiles/tseries.dir/clean:
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/io && $(CMAKE_COMMAND) -P CMakeFiles/tseries.dir/cmake_clean.cmake
.PHONY : lib/io/CMakeFiles/tseries.dir/clean

lib/io/CMakeFiles/tseries.dir/depend:
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/io /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/io /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/io/CMakeFiles/tseries.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/io/CMakeFiles/tseries.dir/depend

