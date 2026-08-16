#pragma once

#include "Acts/Material/MaterialMapper.hpp"
#include "Acts/Utilities/Logger.hpp"

#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Material/IMaterialWriter.hpp"

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
class MaterialMappingAlgorithm : public IAlgorithm {
 public:
  /// @class nested Config class
  /// of the MaterialMapping algorithm
  struct Config {
    /// Input collection
    std::string inputMaterialTracks;
    /// The actually mapped material tracks
    std::string mappedMaterialTracks;
    /// The unmapped part of the material tracks
    std::string unmappedMaterialTracks;
    /// The ACTS material mapper from the core component
    std::shared_ptr<Acts::MaterialMapper> materialMapper;
    /// Material mapping options
    Acts::MaterialMapper::Options materialMapperOptions;
    /// The writer of the material
    std::vector<std::shared_ptr<IMaterialWriter>> materiaMaplWriters;
  };

  /// Constructor
  ///
  /// @param cfg The configuration struct carrying the used tools
  /// @param level The output logging level
  MaterialMappingAlgorithm(const Config& cfg, Acts::Logging::Level level);

  /// Destructor
  /// - it also writes out the file
  ~MaterialMappingAlgorithm() override;

  /// Framework execute method
  ///
  /// @param ctx The algorithm context for event consistency
  ///
  /// @return algorithm process code
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  /// Material mapping state
  std::unique_ptr<Acts::MaterialMapper::State> m_mappingState{nullptr};

  /// Read data handles
  ReadDataHandle<std::vector<Acts::RecordedMaterialTrack>>
      m_inputMaterialTracks{this, "InputMaterialTracks"};

  /// Write data handles
  WriteDataHandle<std::vector<Acts::RecordedMaterialTrack>>
      m_outputMappedMaterialTracks{this, "OutputMappedMaterialTracks"};

  WriteDataHandle<std::vector<Acts::RecordedMaterialTrack>>
      m_outputUnmappedMaterialTracks{this, "OutputUnmappedMaterialTracks"};
};
