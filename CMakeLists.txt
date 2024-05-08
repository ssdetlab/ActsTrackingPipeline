cmake_minimum_required(VERSION 3.14)

# LCG sets CPATH, which gets treated like -I by the compiler. We want to ignore
# warnings from libraries, by unsetting it here, it gets processed by the usual
# target_include_directories call, resulting in the desired -isystem flag.
unset(ENV{CPATH})

project(ActsLUXEPipeline LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

# set the install directories
# to be placed in the build directory
set(CMAKE_INSTALL_BINDIR bin)
set(CMAKE_INSTALL_LIBDIR lib)

# place build products in `<build>/bin` and `<build>/lib` for simple use
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# build options
option(LUXEPipeline_BUILD_EVERYTHING
    "Build with all options enabled" OFF)
option(LUXEPipeline_BUILD_GEOMETRY
    "Build the Geometry component" ON)
option(LUXEPipeline_BUILD_LUXE_READOUT
    "Build the LUXE-specific readout component" OFF)

# TBB, ROOT::Core are required for the
# algorithm execution
find_package(TBB COMPONENTS tbb REQUIRED)
find_package(ROOT COMPONENTS Core REQUIRED)

# Geant4 is required for the geometry
if (LUXEPipeline_BUILD_GEOMETRY)
    find_package(Geant4 REQUIRED)
endif()

# The core components and the Algo
# infrastructure are required
find_package(
    Acts
    COMPONENTS
    Core
    PluginFpeMonitoring
    REQUIRED)

# source files will be added later
add_library(
    ActsLUXEPipelineCore SHARED "")

if (LUXEPipeline_BUILD_EVERYTHING)
    set(LUXEPipeline_BUILD_GEOMETRY ON)
    set(LUXEPipeline_BUILD_LUXE_READOUT ON)
endif()

# collect the library dependencies
add_subdirectory(Infrastructure)
add_subdirectory(TrackFinding)
add_subdirectory(TrackFitting)
add_subdirectory(MagneticField)
add_subdirectory(Material)
add_subdirectory(EventData)
if (LUXEPipeline_BUILD_GEOMETRY)
    add_subdirectory(Geometry)
endif()
if (LUXEPipeline_BUILD_LUXE_READOUT)
    add_subdirectory(Readout)
endif()

# macro to add a run executables
macro(add_run _name _source)
    # automatically prefix the target name
    set(_target "Run${_name}")
    add_executable(${_target} ${_source})
    # core link is public, so no need 
    # to add the includes
    target_link_libraries(
        ${_target}
        PRIVATE
        ActsLUXEPipelineCore)
endmacro()

# run has to be last, so that it can be 
# linked against the library
add_subdirectory(Run)
