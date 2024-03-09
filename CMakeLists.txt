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

find_package(TBB REQUIRED)
find_package(ROOT REQUIRED)
find_package(Geant4 REQUIRED)
find_package(
    Acts
    COMPONENTS
    Core
    PluginGeant4
    PluginFpeMonitoring
    REQUIRED)

message(STATUS "Acts_FOUND: ${Acts_INCLUDE_DIR}")

# source files will be added later
add_library(
    ActsLUXEPipelineCore SHARED "")

# collect the library dependencies
add_subdirectory(Geometry)
#add_subdirectory(EventData)
add_subdirectory(Infrastructure)

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
