target_sources(
    TrackingPipelineCore
    PRIVATE
    RootSimClusterWriter.cpp
    RootSimTrackCandidateWriter.cpp
    RootSimSeedWriter.cpp
    RootFittedSimTrackWriter.cpp
    RootTrackLookupValidationWriter.cpp
    JsonTrackLookupReader.cpp
    JsonTrackLookupWriter.cpp
    RootMaterialTrackReader.cpp  
    RootMaterialTrackWriter.cpp)
