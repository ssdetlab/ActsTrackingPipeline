# 
# Simple run macro is defined upstream
# If needed, add a new run analogously
#

add_run(LookupDataGeneration LookupDataGenerationRun.cpp)
add_run(PathSeeding PathSeedingRun.cpp)
if (LUXEPipeline_BUILD_LUXE_READOUT)
    add_run(TruthTracking TruthTrackingRun.cpp)
endif()
add_run(MaterialMapping MaterialMappingRun.cpp)
add_run(MaterialValidation MaterialValidationRun.cpp)
add_run(FullTracking FullTrackingRun.cpp)
add_run(Test TestRun.cpp)
