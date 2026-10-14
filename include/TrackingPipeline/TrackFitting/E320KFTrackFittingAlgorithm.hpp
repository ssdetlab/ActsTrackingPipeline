#pragma once

#include "Acts/EventData/TrackParameters.hpp"

#include "TrackingPipeline/EventData/E320IndexSeed.hpp"
#include "TrackingPipeline/EventData/E320IndexTrack.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

namespace E320 {

class E320KFTrackFittingAlgorithm : public IAlgorithm {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Input track candidates
    std::string inputTrackCandidates;
    /// Input track parameters
    std::string inputTrackParameters;
    /// Input source links
    std::string inputSourceLinks;
    /// Input track-specific mag field configurations
    std::string inputMagneticFieldParameters;
    /// Output track container
    std::string outputTrackContainer;
    /// Output track container
    std::string outputTracks;
    /// KF fitter
    const KFFitter& fitter;
    /// Maximum number of steps
    std::size_t maxSteps;
    /// KF extensions
    KFFitterExtensions kfExtensions;
    /// Reference surface
    const Acts::Surface* referenceSurface;
    /// Flag to extract the magnetic fields from the 
    /// event magnetic field context
    bool mctxMagneticFields;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  /// @param level logging level
  E320KFTrackFittingAlgorithm(const Config& config, Acts::Logging::Level level);

  /// @brief Destructor
  ~E320KFTrackFittingAlgorithm() override = default;

  /// @brief The execute method
  ///
  /// @param ctx algorithm context of the event
  ///
  /// @return execution process code
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  /// Read data handles
  ReadDataHandle<E320IndexSeeds> m_inputTrackCandidates{
      this, "InputIndexTrackCandidates"};

  ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_inputTrackParameters{this, "InputTrackParameters"};

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  ReadDataHandle<std::vector<std::shared_ptr<MagneticFieldStore>>>
      m_inputMagneticFieldParameters{this, "InputMagneticFieldParameters"};

  /// Write data handles
  WriteDataHandle<KFFitterTrackContainer> m_outputTrackContainer{
      this, "OutputTrackContainer"};

  WriteDataHandle<E320IndexTracks> m_outputTracks{this, "OutputTracks"};
};

}  // namespace E320
