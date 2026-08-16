#pragma once

#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialValidater.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Algorithm performing track propagation and material recording
///
/// The algorithm takes external material tracks (obtained from the material
/// recording procedure), extracts the initial track states (position,
/// momentum), initializes propagation with these parameters and performs the
/// propagation through the geometry, recording material interactions at the
/// material surfaces
class MaterialValidationAlgorithm : public IAlgorithm {
 public:
  /// @brief nested configuration struct
  struct Config {
    // Input material tracks
    std::string inputMaterialTracks;
    // Material validater
    std::shared_ptr<Acts::MaterialValidater> materialValidater;
    /// Output material tracks
    std::string outputMaterialTracks;
  };

  /// Constructor
  ///
  /// @param cfg The configuration struct carrying the used tools
  /// @param level The output logging level
  MaterialValidationAlgorithm(const Config& cfg, Acts::Logging::Level level);

  /// Framework execute method
  ///
  /// @param context The algorithm context for event consistency
  ///
  /// @return algorithm process code
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Configuraion
  Config m_cfg;

  /// Read data handles
  ReadDataHandle<std::vector<Acts::RecordedMaterialTrack>>
      m_inputMaterialTracks{this, "InputMaterialTracks"};

  /// Write data handles
  WriteDataHandle<std::vector<Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "OutputMaterialTracks"};
};
