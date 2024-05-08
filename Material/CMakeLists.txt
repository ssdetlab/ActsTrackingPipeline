#
# Include the material component
# to the pipeline core
#

find_package(ROOT REQUIRED COMPONENTS RIO Hist Tree)

find_package(
    Acts
    COMPONENTS
    PluginJson
    REQUIRED)

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
    src/ROOTMaterialTrackReader.cpp  
    src/ROOTMaterialTrackWriter.cpp  
    src/JsonMaterialWriter.cpp
    src/CoreMaterialMapping.cpp
    src/MaterialValidation.cpp)

target_link_libraries(
    ActsLUXEPipelineCore
    PUBLIC
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Tree
    ActsCore
    ActsPluginJson)
