#pragma once

#include "Acts/Material/MaterialMapper.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Material/IMaterialWriter.hpp"

/// @class CoreMaterialMapping
///
/// @brief Initiates and executes material mapping using the MaterialMapper
/// from the core component of ACTS
///
/// By construction, the material mapping needs inter-event information
/// to build the material maps of accumulated single particle views.
/// However, running it in one single event, puts enormous pressure onto
/// the I/O structure.
///
/// It therefore saves the mapping state/cache as a private member variable
/// and is designed to be executed in a single threaded mode.
class CoreMaterialMapping : public IAlgorithm {
 public:
  /// @class nested Config class
  /// of the MaterialMapping algorithm
  struct Config {
    /// Input collection
    std::string inputMaterialTracks = "material_tracks";

    /// The actually mapped material tracks
    std::string mappedMaterialTracks = "mapped_material_tracks";

    /// The unmapped part of the material tracks
    std::string unmappedMaterialTracks = "unmapped_material_tracks";

    /// The ACTS material mapper from the core component
    std::shared_ptr<Acts::MaterialMapper> materialMapper = nullptr;

    /// The writer of the material
    std::vector<std::shared_ptr<IMaterialWriter>> materiaMaplWriters{};
  };

  /// Constructor
  ///
  /// @param cfg The configuration struct carrying the used tools
  /// @param level The output logging level
  CoreMaterialMapping(const Config& cfg,
                      Acts::Logging::Level level = Acts::Logging::INFO);

  /// Destructor
  /// - it also writes out the file
  ~CoreMaterialMapping() override;

  /// Framework execute method
  ///
  /// @param context The algorithm context for event consistency
  ProcessCode execute(const AlgorithmContext& context) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  std::unique_ptr<Acts::MaterialMapper::State> m_mappingState{nullptr};

  ReadDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_inputMaterialTracks{this, "InputMaterialTracks"};

  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_outputMappedMaterialTracks{this, "OutputMappedMaterialTracks"};

  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_outputUnmappedMaterialTracks{this, "OutputUnmappedMaterialTracks"};
};
