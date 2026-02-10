cmake_minimum_required(VERSION 3.14)

unset(ENV{CPATH})

project(TrackingPipeline LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_INSTALL_BINDIR bin)
set(CMAKE_INSTALL_LIBDIR lib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

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
    HistPainter
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

# EUDAQ configuration (package or local build)
# Allow users to disable EUDAQ package detection if it causes trouble
option(USE_EUDAQ_PACKAGE "Try to use a full EUDAQ CMake package" ON)

if (USE_EUDAQ_PACKAGE)
  find_package(
    eudaq
    2.5.2
    COMPONENTS
      stave
      detector_event
    REQUIRED
  )
  include("${eudaq_DIR}/eudaqStaveConfig.cmake")
else()
  message(STATUS "Skipping EUDAQ package detection, using local from_eudaq build")
  set(EUDAQ_INC_DIR "" CACHE PATH "EUDAQ include directory")
  set(EUDAQ_DICT "" CACHE FILEPATH "EUDAQ detector_event dictionary library")
  if (NOT EXISTS "${EUDAQ_INC_DIR}")
    message(FATAL_ERROR "EUDAQ_INC_DIR does not exist: ${EUDAQ_INC_DIR}")
  endif()
  if (NOT EXISTS "${EUDAQ_DICT}")
    message(FATAL_ERROR "EUDAQ_DICT does not exist: ${EUDAQ_DICT}")
  endif()
endif()


# Dictionary for apollonIo detector_event
add_library(DetEventClusteredDict SHARED)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/TrackingPipeline/Preprocessing
)
ROOT_GENERATE_DICTIONARY(
    G__DetEventClusteredDict 
    ./include/TrackingPipeline/Preprocessing/DetectorEvent.hpp
    MODULE 
    DetEventClusteredDict 
    LINKDEF 
    ./include/TrackingPipeline/Preprocessing/LinkDef.hpp
)
target_link_libraries(
    DetEventClusteredDict
    PRIVATE 
    ROOT::RIO
    ROOT::Physics
    ROOT::GeomPainter
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
    DetEventClusteredDict
)

# Remaining dependencies
find_package(Geant4 REQUIRED)
find_package(Boost REQUIRED)

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

add_library(
    TrackingPipelineCore SHARED ""
)

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
    DetEventClusteredDict
    toml11::toml11
    ActsPreprocessing
)

target_include_directories(
    TrackingPipelineCore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

add_subdirectory(src/Geometry)
add_subdirectory(src/Infrastructure)
add_subdirectory(src/Io)
add_subdirectory(src/MagneticField)
add_subdirectory(src/Material)
add_subdirectory(src/Simulation)
add_subdirectory(src/TrackFinding)
add_subdirectory(src/TrackFitting)
add_subdirectory(src/Alignment)
add_subdirectory(src/Run)
add_subdirectory(src/Analysis)
add_subdirectory(src/Preprocessing)

macro(add_run _name _source)
    set(_target "Run${_name}")
    add_executable(${_target} ${_source})
    target_link_libraries(
        ${_target}
        PRIVATE
        TrackingPipelineCore
    )
endmacro()

add_subdirectory(Run)