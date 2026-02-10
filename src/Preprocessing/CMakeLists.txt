
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

# EUDAQ headers
if (USE_EUDAQ_PACKAGE)
  target_include_directories(
    ActsPreprocessing
    PUBLIC
    "${EUDAQ_INCLUDE_DIRS}" 
  )
else()
  target_include_directories(
    ActsPreprocessing
    PUBLIC
    "${EUDAQ_INC_DIR}"
  )
endif()

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
)

# EUDAQ library / dict
if (USE_EUDAQ_PACKAGE)
  target_link_libraries(
    ActsPreprocessing
    PUBLIC
    eudaq::detector_event       
  )
else()
  target_link_libraries(
    ActsPreprocessing
    PUBLIC
    "${EUDAQ_DICT}"
  )
endif()