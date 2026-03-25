target_sources(
    TrackingPipelineCore
    PRIVATE
    AlignmentAlgorithm.cpp  
    LinearAnnealingScheduler.cpp
    detail/AlignmentStoreBuilders.cpp
    detail/AlignmentStoreUpdaterBuilders.cpp
)
    
