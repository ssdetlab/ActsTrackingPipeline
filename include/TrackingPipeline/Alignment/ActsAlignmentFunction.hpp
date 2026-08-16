#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>

#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"

/// @brief Alignment function calling Acts-implemented alignment machinery
/// for chi2 minimization
class ActsAlignmentFunction : public AlignmentAlgorithm::AlignmentFunction {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// Detector geometry instance
    const Acts::Experimental::Detector* detector;
    /// Magnetic field provider
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
    /// KF extensions
    KFFitterExtensions kfExtensions;
    /// Fitter reference surface
    const Acts::Surface* kfReferenceSurface;
    /// Maximum number of fitter propagation steps
    std::size_t maxKFSteps;
    /// Alignment transform updater
    ActsAlignment::AlignmentTransformUpdater alignmentTransformUpdater;
    /// Alignment transform updater
    ActsAlignment::AlignmentParametersSolver alignmentParametersSolver;
    /// Detector elements to be aligned
    std::vector<Acts::DetectorElementBase*> alignedDetElements;
    /// Cutoff value for average chi2/ndf
    double chi2ONdfCutOff;
    /// Cutoff value for delta of average chi2/ndf within a couple of iterations
    std::pair<std::size_t, double> deltaChi2ONdfCutOff;
    /// Maximum number of iterations
    std::size_t maxAlignmentFitNumIt;
    /// Alignment mask
    ActsAlignment::AlignmentMask alignmentMask;
    /// Number of re-fitting iterations
    std::size_t nRefittingIt;
    /// Track parameters estimator
    std::shared_ptr<ITrackParametersEstimator> trackParametersEstimator;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit ActsAlignmentFunction(const Config& cfg);

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
  AlignmentAlgorithm::AlignmentResult operator()(
      const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
      const Acts::CalibrationContext& cctx,
      const AlignmentAlgorithm::SourceLinkContainer& alignmentFitSourceLinks,
      const AlignmentAlgorithm::SourceLinkContainer&
          initialTrackStateFitSourceLinks,
      const AlignmentAlgorithm::TrackParametersContainer& initialParameters,
      const AlignmentAlgorithm::MagneticFieldParametersContainer&
          magFieldParameters) const override;

 private:
  /// Configuration
  Config m_cfg;

  /// Acts alignment instance
  std::unique_ptr<ActsAligner> m_align;
};
