#
# Include the infrastructure component
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
    src/IAlgorithm.cpp
    src/Sequencer.cpp
    src/WhiteBoard.cpp
    src/SequenceElement.cpp)

# For now we assume that the pipeline 
# always uses the Geant4 plugin
target_link_libraries(
    ActsLUXEPipelineCore
    PUBLIC
    TBB::tbb
    ROOT::Core
    ActsCore
    ActsPluginFpeMonitoring)
