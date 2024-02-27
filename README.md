# ActsLUXEPipeline

-- Description
Implementation of the track reconstruction pipeline for the LUXE experiment.

The modules are compiled in a single shared library to be linked 
against the executables. 

The build produces a set of executable Runs that can be found in the 
bin folder of the build directory. 

The available pipelines are:

-- EnergyPositionHistosRun
Create a set of histograms and fits to be used in the seeding

-- Installation guide
To install the library one needs to have Acts Core and Geant4 plugin
compiled and installed locally. DON'T FORGET TO RUN "this_acts.sh" to 
ensure the environment variables are in place. The build will fail otherwise. 
Otherwise, the build process follows the standard scheme. 