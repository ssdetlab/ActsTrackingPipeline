target_sources(
    TrackingPipelineCore
    PRIVATE
    RootSimClusterWriter.cpp
    RootSimSeedWriter.cpp
    RootSimTrackCandidateWriter.cpp
    RootSimTrackWriter.cpp
    RootMaterialTrackReader.cpp  
    RootMaterialTrackWriter.cpp
    RootTrackLookupValidationWriter.cpp
    RootTrackParamsReader.cpp
    RootSimDataReader.cpp
    JsonTrackLookupReader.cpp
    JsonTrackLookupWriter.cpp
)
