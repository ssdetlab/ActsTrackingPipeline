#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "TrackingPipeline/Alignment/IAnnealingScheduler.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"
#include "TrackingPipeline/Utilities/ConcatVectorView.hpp"
#include "TrackingPipeline/Utilities/NonOwningProjectionContainer.hpp"

class AlignmentAlgorithm final : public IAlgorithm {
 public:
  using SourceLinkContainer = std::vector<
      NonOwningProjectionContainer<ConcatVectorView<Acts::SourceLink>>>;

  using TrackParametersContainer = NonOwningProjectionContainer<
      std::vector<Acts::CurvilinearTrackParameters>>;

  /// Alignment function that takes the above parameters and runs alignment
  class AlignmentFunction {
   public:
    virtual ~AlignmentFunction() = default;
    virtual AlignmentResult operator()(
        const SourceLinkContainer& sourceLinks,
        const TrackParametersContainer& initialParameters,
        const ActsAlignment::AlignmentOptions<KFFitterOptions>& options) = 0;
  };

  /// Create the alignment function implementation.
  static std::shared_ptr<AlignmentFunction> makeAlignmentFunction(
      const std::shared_ptr<const Acts::Experimental::Detector>& detector,
      const std::shared_ptr<const Acts::MagneticFieldProvider>& magneticField);

  /// Config struct
  struct Config {
    /// Input source links
    std::string inputSourceLinks;
    /// Input track candidates
    std::string inputTrackCandidates;
    /// Input track parameters
    std::string inputTrackParameters;
    /// Output alignment parameters
    std::string outputAlignmentParameters;
    /// Output updated track parameters
    std::string outputTrackParameters;
    /// Type erased fitter function.
    std::shared_ptr<AlignmentFunction> alignmentFunction;
    /// Alignment transform updater
    ActsAlignment::AlignmentTransformUpdater alignmentTransformUpdater;
    /// Alignment transform updater
    ActsAlignment::AlignmentParametersSolver alignmentParametersSolver;
    /// Detector elements to be aligned
    std::vector<Acts::DetectorElementBase*> alignedDetElements;
    /// KF options
    KFFitterOptions kfOptions;
    /// Cutoff value for average chi2/ndf
    double chi2ONdfCutOff;
    /// Cutoff value for delta of average chi2/ndf within a couple of iterations
    std::pair<std::size_t, double> deltaChi2ONdfCutOff;
    /// Maximum number of iterations
    std::size_t maxAlignmentFitNumIt;
    /// Alignment mask
    ActsAlignment::AlignmentMask alignmentMask;
    /// Initial track state covariance prior
    Acts::BoundMatrix originCov;
    /// Constraints vector
    std::vector<Acts::SourceLink> constraints;
    /// Number of iterations of the fit with initial guess refitting
    std::size_t nRefittingIt;
    /// Anneling factor scheduler
    std::shared_ptr<IAnnealingScheduler> annealingScheduler;
    /// Track parameters estimator
    std::shared_ptr<ITrackParametersEstimator> trackParametersEstimator;
  };

  /// Constructor of the alignment algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  AlignmentAlgorithm(const Config& cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the alignment algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ProcessCode execute(const AlgorithmContext& ctx) const override;

 private:
  Config m_cfg;

  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  ReadDataHandle<IndexSeeds> m_inputTrackCandidates{this,
                                                    "InputTrackCandidates"};

  ReadDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_inputTrackParameters{this, "InputTrackParameters"};

  WriteDataHandle<AlignmentResult> m_outputAlignmentParameters{
      this, "OutputAlignmentParameters"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParameters{this, "OutputTrackParameters"};
};
