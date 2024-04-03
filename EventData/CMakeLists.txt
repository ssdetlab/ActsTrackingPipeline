#
# Include the event data component
# to the pipeline core
#

# ROOT components for the 
# dictionary generation
# and event data reading

#find_package(
#    ROOT
#    COMPONENTS
#    RIO
#    Tree
#    Physics
#    REQUIRED)

# Includes for Acts are propagated
# through the library linking
target_include_directories(
    ActsLUXEPipelineCore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(
    ActsLUXEPipelineCore
    PUBLIC
    ActsCore)
      