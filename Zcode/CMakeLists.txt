#
# Include the Utilities components
# to the pipeline core
#

find_package(ROOT COMPONENTS Tree Hist REQUIRED)

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
    ROOT::Core
    ROOT::Imt
    ROOT::RIO
    ROOT::Net
    ROOT::Hist
    ROOT::Graf
    ROOT::Graf3d
    ROOT::Gpad
    ROOT::ROOTVecOps
    ROOT::Tree
    ROOT::TreePlayer
    ROOT::Rint
    ROOT::Postscript
    ROOT::Matrix
    ROOT::Physics
    ROOT::MathCore
    ROOT::Thread
    ROOT::MultiProc
    ROOT::MultiProc
    ROOT::ROOTDataFrame
)
      