cmake_minimum_required(VERSION 3.14)

# LCG sets CPATH, which gets treated like -I by the compiler. We want to ignore
# warnings from libraries, by unsetting it here, it gets processed by the usual
# target_include_directories call, resulting in the desired -isystem flag.
unset(ENV{CPATH})

project(TrackingPipeline LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# set the install directories
# to be placed in the build directory
set(CMAKE_INSTALL_BINDIR bin)
set(CMAKE_INSTALL_LIBDIR lib)

# place build products in `<build>/bin` and `<build>/lib` for simple use
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# TBB, ROOT::Core are required for the
# algorithm execution
find_package(
    TBB 
    COMPONENTS 
    tbb 
    REQUIRED
)
find_package(
    ROOT 
    COMPONENTS 
    Core
    RIO
    Tree
    Hist
    Gpad
    Net
    Physics
    GeomPainter
    REQUIRED
)

find_package(
  OpenMP 
  COMPONENTS
  CXX
  REQUIRED
)

find_package(
  eudaq_preprocessor
  COMPONENTS
  detector_event
  REQUIRED
)

# Dictionary for the sim data
add_library(SimEventDict SHARED)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/TrackingPipeline/Io/detail
)
ROOT_GENERATE_DICTIONARY(
    G__SimEventDict
    ./include/TrackingPipeline/Io/detail/ROOTDictDefs.h
    MODULE 
    SimEventDict 
    LINKDEF 
    ./include/TrackingPipeline/Io/detail/LinkDef.h
)
target_link_libraries(
    SimEventDict 
    PRIVATE 
    ROOT::RIO
    ROOT::Physics
    ROOT::GeomPainter
    eudaq_preprocessor::detector_event
)

# Geant4 is required for the geometry
find_package(Geant4 REQUIRED)

# Boost is required for the Acts
find_package(Boost REQUIRED)

# The core components and the Algo
# infrastructure are required
find_package(
    Acts
    COMPONENTS
    Core
    Fatras
    PluginGeant4
    PluginJson
    PluginFpeMonitoring
    Alignment
    REQUIRED
)

find_package(toml11 REQUIRED)

# source files will be added later
add_library(
    TrackingPipelineCore SHARED "")

target_link_libraries(
    TrackingPipelineCore
    PUBLIC
    TBB::tbb
    OpenMP::OpenMP_CXX
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Physics
    ROOT::Hist
    ROOT::Gpad
    ROOT::Net
    G__SimEventDict
    ActsCore
    ActsFatras
    ActsPluginGeant4
    ActsPluginJson
    ActsPluginFpeMonitoring
    ActsAlignment
    eudaq_preprocessor::detector_event
    toml11::toml11
)

target_include_directories(
    TrackingPipelineCore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# collect the library dependencies
add_subdirectory(src/Geometry)
add_subdirectory(src/Infrastructure)
add_subdirectory(src/Io)
add_subdirectory(src/MagneticField)
add_subdirectory(src/Material)
add_subdirectory(src/Simulation)
add_subdirectory(src/TrackFinding)
add_subdirectory(src/TrackFitting)
add_subdirectory(src/Alignment)

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
        TrackingPipelineCore
    )
endmacro()

# run has to be last, so that it can be 
# linked against the library
add_subdirectory(Run)
