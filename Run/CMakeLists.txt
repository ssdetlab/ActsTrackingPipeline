# 
# Simple run macro is defined upstream
# If needed, add a new run analogously
#

add_run(LookupDataGeneration LookupDataGenerationRun.cpp)
add_run(LookupDataValidation LookupDataValidationRun.cpp)
add_run(MaterialMapping MaterialMappingRun.cpp)
add_run(MaterialValidation MaterialValidationRun.cpp)
add_run(FullTracking FullTrackingRun.cpp)
add_run(FastSim FastSimRun.cpp)
