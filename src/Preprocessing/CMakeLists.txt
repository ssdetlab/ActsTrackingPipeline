
add_library(
    ActsPreprocessing SHARED
    AnalysisFunctions.cpp
    # Preprocessing.cpp
)

target_include_directories(
    ActsPreprocessing
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
)

target_link_libraries(
    ActsPreprocessing
    PUBLIC
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Gpad
    ROOT::Physics
    DetEventClusteredDict
    DetEventDict
)