#pragma once

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"

/// @brief alignment function perform global Minuit-based optimization
///
/// Function passed to the alignment algorithm performing transverse shift
/// and rotation about the COM horizontal axis alignment search using
/// MINUIT chi2 minimization
class MinuitGlobalAlignmentFunction
    : public AlignmentAlgorithm::AlignmentFunction {
 public:
  /// @brief nested configuration struct
  struct Config {
    /// Detector geometry instance
    const Acts::Experimental::Detector* detector;
    /// Magnetic field provider
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
    /// Detector elements to be aligned
    std::vector<Acts::DetectorElementBase*> alignedDetElements;
    /// Maximum number of fitter propagation steps
    std::size_t maxKFSteps;
    /// KF extensions
    KFFitterExtensions kfExtensions;
    /// Fitter reference surface
    const Acts::Surface* kfReferenceSurface;
    /// Track parameters estimator
    std::shared_ptr<ITrackParametersEstimator> trackParametersEstimator;
  };

  /// @brief MINUIT FCN implementation
  class MinuitAlignmentFunctionImpl : public ROOT::Minuit2::FCNBase {
   public:
    /// @brief nested options struct
    struct Options {
      /// Source links used in the track fit optimization
      const AlignmentAlgorithm::SourceLinkContainer&
          trackFitSourceLinkContainer;
      /// Source links used in the inital track state estimation
      const AlignmentAlgorithm::SourceLinkContainer&
          trackParametersFitSourceLinkContainer;
      /// Initial track parameters estimates
      const AlignmentAlgorithm::TrackParametersContainer& trackParameters;
      /// Track-specific magnetic fields configuration
      const AlignmentAlgorithm::MagneticFieldParametersContainer&
          magFieldParameters;
      /// Kalman filter instance
      const KFFitter& fitter;
      /// Kalman filter extensions
      const KFFitterExtensions& kfExtensions;
      /// Kalman filter reference surface
      const Acts::Surface* kfReferenceSurface;
      /// Maximum number of KF steps
      std::size_t maxKFSteps;
      /// Initial track state estimator
      std::shared_ptr<ITrackParametersEstimator> trackParametersEstimator;
      /// Detector geometry instance
      const Acts::Experimental::Detector* detector;
    };

    explicit MinuitAlignmentFunctionImpl(const Options& opt);

    double Up() const override;

    double operator()(const std::vector<double>& pars) const override;

   private:
    Options m_opt;
  };

  explicit MinuitGlobalAlignmentFunction(const Config& cfg);

  AlignmentResult operator()(
      const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
      const Acts::CalibrationContext& cctx,
      const AlignmentAlgorithm::SourceLinkContainer& sourceLinks,
      const AlignmentAlgorithm::TrackParametersContainer& initialParameters,
      const AlignmentAlgorithm::MagneticFieldParametersContainer&
          magFieldParameters) override;

 private:
  Config m_cfg;
};
