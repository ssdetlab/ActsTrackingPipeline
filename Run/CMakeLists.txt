# 
# Simple run macro is defined upstream
# If needed, add a new run analogously
#

add_run(EnergyPositionHistos EnergyPositionHistosRun.cpp)
if (LUXEPipeline_BUILD_LUXE_READOUT)
    add_run(TruthTracking TruthTrackingRun.cpp)
endif()
