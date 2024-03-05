#
# Include the event data component
# to the pipeline core
#

# Includes for Acts are propagated
# through the library linking
target_include_directories(
    ActsLUXEPipelineCore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

# Sources
target_sources(
    ActsLUXEPipelineCore
    PRIVATE
    src/ROOTDataReader.cpp)

# Dictionary for the sim data
add_library(LUXESimEventDict SHARED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
ROOT_GENERATE_DICTIONARY(
    G__LUXESimEventDict 
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ActsLUXEPipeline/ROOTDictDefs.h
    MODULE 
    LUXESimEventDict 
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/ActsLUXEPipeline/LinkDef.h)
target_link_libraries(
    G__LUXESimEventDict 
    PRIVATE 
    ROOT::RIO
    ROOT::Physics)

# For now we assume that the pipeline 
# always uses the Geant4 plugin
target_link_libraries(
    ActsLUXEPipelineCore
    PUBLIC
    ActsCore
    ROOT::Tree
    ROOT::Physics
    G__LUXESimEventDict)
      