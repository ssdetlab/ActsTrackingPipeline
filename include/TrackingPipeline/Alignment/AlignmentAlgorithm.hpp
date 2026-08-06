#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <memory>
#include <unordered_set>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"
#include "TrackingPipeline/Utilities/NonOwningProjectionContainer.hpp"

/// @brief Algorihtm perfroming tracking detector local and global alignment
class AlignmentAlgorithm final : public IAlgorithm {
 public:
  /// Containers shorthands
  using SourceLinkContainer = std::vector<
      detail::NonOwningProjectionContainer<std::vector<Acts::SourceLink>>>;
  using TrackParametersContainer = detail::NonOwningProjectionContainer<
      std::vector<Acts::CurvilinearTrackParameters>>;
  using MagneticFieldParametersContainer = detail::NonOwningProjectionContainer<
      std::vector<std::shared_ptr<MagneticFieldStore>>>;

  /// @brief Alignment function wrapper that takes the above parameters
  /// and runs the alignment procedure
  class AlignmentFunction {
   public:
    /// @brief operator performing the alignment procedure
    ///
    /// @param gctx current geometry context
    /// @param mctx current magnetic field context
    /// @param cctx current calibration context
    /// @param trackFitSourceLinks source links participating in the alignment fit
    /// @param initialTrackStateFitSourceLinks source links participating in
    /// the initial track state fit
    /// @param initialParameters initial track parameters
    /// @param magFieldParameters track-specific magnetic field configurations
    ///
    /// @return struct containing the results of the alignment procedure
    virtual AlignmentResult operator()(
        const Acts::GeometryContext& gctx,
        const Acts::MagneticFieldContext& mctx,
        const Acts::CalibrationContext& cctx,
        const AlignmentAlgorithm::SourceLinkContainer& alignmentFitSourceLinks,
        const AlignmentAlgorithm::SourceLinkContainer&
            initialTrackStateFitSourceLinks,
        const AlignmentAlgorithm::TrackParametersContainer& initialParameters,
        const AlignmentAlgorithm::MagneticFieldParametersContainer&
            magFieldParameters) = 0;
  };

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
  WriteDataHandle<AlignmentResult> m_outputAlignmentParameters{
      this, "OutputAlignmentParameters"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParameters{this, "OutputTrackParameters"};
};
