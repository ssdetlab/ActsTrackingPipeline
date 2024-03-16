#
# Include the readout component
# to the pipeline core
#

# ROOT components for the 
# dictionary generation
# and event data reading
find_package(
    ROOT 
    COMPONENTS 
    RIO 
    Tree
    Physics
    REQUIRED)

# Includes for Acts are propagated
# through the library linking
target_include_directories(
    ActsLUXEPipelineCore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Sources
# Note: No sources for now
# target_sources(
    # ActsLUXEPipelineCore
    # PRIVATE)

# Dictionary for the sim data
add_library(LUXESimEventDict SHARED)
message(STATUS "ROOT_INCLUDE_DIRS: ${CMAKE_CURRENT_SOURCE_DIR}")
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ActsLUXEPipeline/detail)
ROOT_GENERATE_DICTIONARY(
    G__LUXESimEventDict
    ./include/ActsLUXEPipeline/detail/ROOTDictDefs.h
    MODULE 
    LUXESimEventDict 
    LINKDEF 
    ./include/ActsLUXEPipeline/detail/LinkDef.h)
target_link_libraries(
    G__LUXESimEventDict 
    PRIVATE 
    ROOT::RIO
    ROOT::Physics)

target_link_libraries(
    ActsLUXEPipelineCore
    PUBLIC
    ActsCore
    ROOT::Tree
    ROOT::Physics
    G__LUXESimEventDict)
      