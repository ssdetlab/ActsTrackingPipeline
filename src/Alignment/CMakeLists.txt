target_sources(
    TrackingPipelineCore
    PRIVATE
    AlignmentAlgorithm.cpp  
    LinearAnnealingScheduler.cpp
    LocalAlignmentParametersSolverConstraints.cpp
    LocalAlignmentParametersSolverSVD.cpp
    GlobalAlignmentParametersSolver.cpp
    detail/AlignmentStoreBuilders.cpp
    detail/AlignmentStoreUpdaterBuilders.cpp
)
    
