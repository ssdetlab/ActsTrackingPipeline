#
# Include the Utilities components
# to the pipeline core
#

# FOR TESTING
find_package(Boost REQUIRED COMPONENTS serialization)

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
)

target_link_libraries(
    ActsLUXEPipelineCore
    PUBLIC
    ActsCore
)
      