#pragma once

#include <string>
#include <vector>

namespace TrackingPipeline::Preprocessing {

struct PreprocessingConfig {
  // Directories containing input ROOT files (EUDAQ ev_data folders)
  std::vector<std::string> inputDirs;

  // Name of the input tree inside those files (e.g. "MyTree")
  std::string inputTreeName;

  // Output ROOT file with ApollonIo::DetectorEvent tree
  std::string outputFile;

  // Branch name for the input detector_event object (usually "event")
  std::string inputBranchName = "event";

  // Name of the output tree (e.g. "MyTree")
  std::string outputTreeName = "MyTree";

  // Number of initial entries to skip (for debugging or warm‑up)
  std::size_t skipEntries = 0;
};

// Run the preprocessing: read EUDAQ detector_event, clean & cluster, write ApollonIo::DetectorEvent.
void runPreprocessing(const PreprocessingConfig& cfg);

}  // namespace TrackingPipeline::Preprocessing
