#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <mutex>
#include <string>
#include <vector>

#include "TChain.h"
#include "TrackingPipeline/Io/ITrackParamsReader.hpp"

class RootIPPositronTrackParametersReader final : public ITrackParamsReader {
 public:
  struct Config {
    /// Input file path
    std::vector<std::string> filePaths;
    /// Input tree name
    std::string treeName;
    Acts::Transform3 transform;
  };

  RootIPPositronTrackParametersReader(const Config& config);

  std::vector<Acts::CurvilinearTrackParameters> read() override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  TChain* m_chain = nullptr;

  std::vector<float>* m_positionX = nullptr;
  std::vector<float>* m_positionY = nullptr;
  std::vector<float>* m_positionZ = nullptr;

  std::vector<float>* m_phi = nullptr;
  std::vector<float>* m_theta = nullptr;
  std::vector<float>* m_E = nullptr;

  std::vector<int>* m_pdgId = nullptr;

  std::mutex m_readMutex;
};
