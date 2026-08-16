#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMigrad.h"
#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"
#include "TrackingPipeline/Alignment/AlignmentStore.hpp"
#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"
#include "TrackingPipeline/TrackFinding/ITrackParametersEstimator.hpp"

namespace E320 {

/// @brief e320-specific alignment function performing global Minuit-based optimization
///
/// Function passed to the alignment algorithm performing transverse shift
/// and rotation about the COM horizontal axis alignment search using
/// MINUIT chi2 minimization
class E320MinuitGlobalAlignmentFunction
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
    /// Alignment store instance
    std::shared_ptr<AlignmentStore> alignmentStore;
    /// KF extensions
    KFFitterExtensions kfExtensions;
    /// Fitter reference surface
    const Acts::Surface* kfReferenceSurface;
    /// Maximum number of fitter propagation steps
    std::size_t maxKFSteps;
    /// Track parameters estimator
    std::shared_ptr<ITrackParametersEstimator> trackParametersEstimator;
    /// Initial parameters in global coordinates [dy, dz, dtheta]
    Acts::ActsDynamicVector initialPars;
    /// Initial parameters' steps
    Acts::ActsDynamicVector initialParSteps;
    /// Parameters' low bounds
    Acts::ActsDynamicVector parLowBounds;
    /// Parameters' high bounds
    Acts::ActsDynamicVector parHighBounds;
    /// Up level of the FNC
    double upLevel;
  };

  /// @brief MINUIT FCN implementation
  class MinuitAlignmentFunctionImpl : public ROOT::Minuit2::FCNBase {
   public:
    /// @brief nested options struct
    struct Options {
      /// Geometry context
      const Acts::GeometryContext& gctx;
      /// Magnetic field context
      const Acts::MagneticFieldContext& mctx;
      /// Calibration context
      const Acts::CalibrationContext& cctx;
      /// Source links used in the track fit optimization
      const AlignmentAlgorithm::SourceLinkContainer&
          alignmentFitSourceLinkContainer;
      /// Source links used in the inital track state estimation
      const AlignmentAlgorithm::SourceLinkContainer&
          initialTrackStateFitSourceLinkContainer;
      /// Initial track parameters estimates
      const AlignmentAlgorithm::TrackParametersContainer& trackParameters;
      /// Track-specific magnetic fields configuration
      const AlignmentAlgorithm::MagneticFieldParametersContainer&
          magFieldParameters;
      /// Kalman filter instance
      const std::shared_ptr<KFFitter>& fitter;
      /// Function config instance
      const Config& functionCfg;
    };

    /// @brief nested state struct
    struct State {
      /// Current chi2 estimate
      double chi2;
      /// Number of the measurement d.o.fs
      std::size_t measurementDim;
    };

    /// @brief Constuctor
    ///
    /// @param opt options struct
    /// @param level logging level
    MinuitAlignmentFunctionImpl(const Options& opt, Acts::Logging::Level level);

    /// @brief get const access to the function state
    const State& state() const;

    /// @brief return FCN up level
    double Up() const override;

    /// @brief perform parameter optimization
    ///
    /// @param pars vector of MINUIT parameters
    ///
    /// @return chi2 post optimization step
    double operator()(const std::vector<double>& pars) const override;

   protected:
    /// @brief access to logging instance
    const Acts::Logger& logger() const { return *m_logger; }

   private:
    /// Options
    Options m_opt;

    /// Current function state
    std::unique_ptr<State> m_state;

    /// Logging instance
    std::unique_ptr<const Acts::Logger> m_logger;

    /// Alignment surfaces COM
    Acts::Vector3 m_com;

    /// Alignment surfaces lever arms
    std::unordered_map<Acts::GeometryIdentifier, Acts::Vector3> m_leverArms;

    /// Nominal geometry context
    Acts::GeometryContext m_nominalGctx;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  /// @param level logging level
  E320MinuitGlobalAlignmentFunction(const Config& cfg,
                                    Acts::Logging::Level level);

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

 protected:
  /// @brief access to logging instance
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  /// Configuration
  Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Logging level
  Acts::Logging::Level m_level;

  /// Kalman filter instance
  std::shared_ptr<KFFitter> m_fitter;
};

}  // namespace E320
