#pragma once

#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialValidater.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @class MaterialValidation
///
class MaterialValidation : public IAlgorithm {
 public:
  /// @class nested Config class
  /// of the MaterialMapping algorithm
  struct Config {
    // Input file with validation tracks
    std::string inputMaterialTracks = "inputMaterialTracks";
    // The validater
    std::shared_ptr<Acts::MaterialValidater> materialValidater = nullptr;
    /// Output collection name
    std::string outputMaterialTracks = "material_tracks";
  };

  /// Constructor
  ///
  /// @param cfg The configuration struct carrying the used tools
  /// @param level The output logging level
  MaterialValidation(const Config& cfg,
                     Acts::Logging::Level level = Acts::Logging::INFO);

  /// Destructor
  /// - it also writes out the file
  ~MaterialValidation() override = default;

  /// Framework execute method
  ///
  /// @param context The algorithm context for event consistency
  ProcessCode execute(const AlgorithmContext& context) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Configuraion
  Config m_cfg;

  ReadDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_inputMaterialTracks{this, "InputMaterialTracks"};

  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "OutputMaterialTracks"};
};
