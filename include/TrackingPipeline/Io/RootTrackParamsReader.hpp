#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <mutex>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TrackingPipeline/Io/ITrackParamsReader.hpp"

class RootTrackParamsReader final : public ITrackParamsReader {
 public:
  struct TrackParams {
    double positionX;
    double positionY;
    double positionZ;

    double phi;
    double theta;

    double qOverP;
    int pdgId;
  };

  struct Config {
    /// Input file path
    std::vector<std::string> filePaths;
    /// Input tree name
    std::string treeName;
    Acts::Transform3 transform;
  };

  explicit RootTrackParamsReader(const Config& config);

  std::vector<Acts::CurvilinearTrackParameters> read() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  TFile* m_file = nullptr;
  TTree* m_tree = nullptr;
  TChain* m_chainOwner = nullptr;

  TrackParams m_params{};

  std::mutex m_readMutex;
};
