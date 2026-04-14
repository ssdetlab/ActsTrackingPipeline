# 
# Simple run macro is defined upstream
# If needed, add a new run analogously
#

add_run(MaterialMapping MaterialMappingRun.cpp)
add_run(MaterialValidation MaterialValidationRun.cpp)
# add_run(LookupDataGeneration LookupDataGenerationRun.cpp)
# add_run(LookupDataValidation LookupDataValidationRun.cpp)
add_run(FullTracking FullTrackingRun.cpp)
add_run(FullBpmTracking FullBpmTrackingRun.cpp)
add_run(AlignmentLocal AlignmentLocalRun.cpp)
add_run(AlignmentGlobal AlignmentGlobalRun.cpp)
add_run(FastSim FastSimRun.cpp)
add_run(FullSimTracking FullSimTrackingRun.cpp)
add_run(FullSimBpmTracking FullSimBpmTrackingRun.cpp)
add_run(SimAlignmentLocal SimAlignmentLocalRun.cpp)
add_run(SimAlignmentGlobal SimAlignmentGlobalRun.cpp)
