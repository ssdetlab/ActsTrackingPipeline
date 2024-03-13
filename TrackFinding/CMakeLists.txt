#
# Include the track finding component
# to the pipeline core
#

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

target_link_libraries(
    ActsLUXEPipelineCore
    PUBLIC
    ActsCore)
      