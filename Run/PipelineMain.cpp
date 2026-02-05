#include "TrackingPipeline/Run/PipelineRun.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <Pipeline.conf>\n";
    return 1;
  }
  return TrackingPipeline::runPipeline(argv[1]);
}
