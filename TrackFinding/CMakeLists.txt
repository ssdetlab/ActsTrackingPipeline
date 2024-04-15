find_package(ROOT COMPONENTS Tree REQUIRED)

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
        ROOT::Tree
)
      