target_sources(
    TrackingPipelineCore
    PRIVATE
    AlignmentAlgorithm.cpp  
    LinearAnnealingScheduler.cpp
    LocalAlignmentParametersSolverConstraints.cpp
    LocalAlignmentParametersSolverSVD.cpp
    GlobalAlignmentParametersSolver.cpp
    LocalAlignmentTransformUpdater.cpp
    GlobalAlignmentTransformUpdater.cpp
    detail/AlignmentStoreBuilders.cpp
)
    
