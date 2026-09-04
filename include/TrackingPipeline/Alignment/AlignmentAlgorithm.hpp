#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include <memory>
#include <unordered_set>
#include <vector>

#include "TrackingPipeline/Alignment/AlignmentFunction.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

/// @brief Algorihtm perfroming tracking detector local and global alignment
///
/// Algorithm takes in seeds with estimated initial track states and performs
/// alignment of the tracking detector. The alignment procedures should be
/// implemented in the children of the type-erased AlignmentFunction nested
/// class. It is assumed that the alignment function takes in contexts and
/// colletions of source links, track parameters, and magnetic fields required
/// for track fitting and initial track state re-esimation.
class AlignmentAlgorithm final : public IAlgorithm {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// Input source links
    std::string inputSourceLinks;
    /// Input track candidates
    std::string inputTrackCandidates;
    /// Input track parameters
    std::string inputTrackParameters;
    /// Input track-specific mag field configurations
    std::string inputMagneticFieldParameters;
    /// Output alignment parameters
    std::string outputAlignmentParameters;
    /// Output updated track parameters
    std::string outputTrackParameters;
    /// Type erased alignment fit function
    std::shared_ptr<AlignmentFunction> alignmentFunction;
    /// Range of surfaces used for the track fit
    std::unordered_set<Acts::GeometryIdentifier> alignmentFitSurfaces;
    /// Range of surfaces used for the initial track state fit
    std::unordered_set<Acts::GeometryIdentifier> initialTrackStateFitSurfaces;
  };

  /// @brief Constructor
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  AlignmentAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief execute method
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  ///
  /// @return a process code to steer the algorithm flow
  ProcessCode execute(const AlgorithmContext& ctx) const override;

 private:
  /// Configuration
  Config m_cfg;

  /// Read data handles
  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  ReadDataHandle<IndexSeeds> m_inputTrackCandidates{this,
                                                    "InputTrackCandidates"};

  ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_inputTrackParameters{this, "InputTrackParameters"};

  ReadDataHandle<std::vector<std::shared_ptr<MagneticFieldStore>>>
      m_inputMagneticFieldParameters{this, "InputMagneticFieldParameters"};

  /// Write data handles
  WriteDataHandle<AlignmentFunction::AlignmentResult>
      m_outputAlignmentParameters{this, "OutputAlignmentParameters"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParameters{this, "OutputTrackParameters"};
};
