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
include lib/parallel/CMakeFiles/parallel.dir/depend.make

# Include the progress variables for this target.
include lib/parallel/CMakeFiles/parallel.dir/progress.make

# Include the compile flags for this target's objects.
include lib/parallel/CMakeFiles/parallel.dir/flags.make

lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o: lib/parallel/CMakeFiles/parallel.dir/flags.make
lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o: ../../lib/parallel/parallel.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel.dir/parallel.cc.o -c /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/parallel/parallel.cc

lib/parallel/CMakeFiles/parallel.dir/parallel.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel.dir/parallel.cc.i"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/parallel/parallel.cc > CMakeFiles/parallel.dir/parallel.cc.i

lib/parallel/CMakeFiles/parallel.dir/parallel.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel.dir/parallel.cc.s"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/parallel/parallel.cc -o CMakeFiles/parallel.dir/parallel.cc.s

lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o.requires:

.PHONY : lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o.requires

lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o.provides: lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o.requires
	$(MAKE) -f lib/parallel/CMakeFiles/parallel.dir/build.make lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o.provides.build
.PHONY : lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o.provides

lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o.provides.build: lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o


lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o: lib/parallel/CMakeFiles/parallel.dir/flags.make
lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o: ../../lib/parallel/mpidata.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel.dir/mpidata.cc.o -c /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/parallel/mpidata.cc

lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel.dir/mpidata.cc.i"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/parallel/mpidata.cc > CMakeFiles/parallel.dir/mpidata.cc.i

lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel.dir/mpidata.cc.s"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && /opt/software/intel_2015.u2/impi/5.0.3.048/intel64/bin/mpicxx  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/parallel/mpidata.cc -o CMakeFiles/parallel.dir/mpidata.cc.s

lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o.requires:

.PHONY : lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o.requires

lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o.provides: lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o.requires
	$(MAKE) -f lib/parallel/CMakeFiles/parallel.dir/build.make lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o.provides.build
.PHONY : lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o.provides

lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o.provides.build: lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o


# Object files for target parallel
parallel_OBJECTS = \
"CMakeFiles/parallel.dir/parallel.cc.o" \
"CMakeFiles/parallel.dir/mpidata.cc.o"

# External object files for target parallel
parallel_EXTERNAL_OBJECTS =

lib/parallel/libparallel.a: lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o
lib/parallel/libparallel.a: lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o
lib/parallel/libparallel.a: lib/parallel/CMakeFiles/parallel.dir/build.make
lib/parallel/libparallel.a: lib/parallel/CMakeFiles/parallel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libparallel.a"
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && $(CMAKE_COMMAND) -P CMakeFiles/parallel.dir/cmake_clean_target.cmake
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/parallel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/parallel/CMakeFiles/parallel.dir/build: lib/parallel/libparallel.a

.PHONY : lib/parallel/CMakeFiles/parallel.dir/build

lib/parallel/CMakeFiles/parallel.dir/requires: lib/parallel/CMakeFiles/parallel.dir/parallel.cc.o.requires
lib/parallel/CMakeFiles/parallel.dir/requires: lib/parallel/CMakeFiles/parallel.dir/mpidata.cc.o.requires

.PHONY : lib/parallel/CMakeFiles/parallel.dir/requires

lib/parallel/CMakeFiles/parallel.dir/clean:
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel && $(CMAKE_COMMAND) -P CMakeFiles/parallel.dir/cmake_clean.cmake
.PHONY : lib/parallel/CMakeFiles/parallel.dir/clean

lib/parallel/CMakeFiles/parallel.dir/depend:
	cd /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/lib/parallel /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel /home/anas/Anas/saras/latest_saras/Omega_per_g/time_average/new/compile/build/lib/parallel/CMakeFiles/parallel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/parallel/CMakeFiles/parallel.dir/depend
