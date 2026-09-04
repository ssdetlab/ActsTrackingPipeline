#include "TrackingPipeline/Alignment/ActsAlignmentFunction.hpp"

#include <cstddef>

#include "TrackingPipeline/Infrastructure/TypeDefinitions.hpp"

ActsAlignmentFunction::ActsAlignmentFunction(const Config& cfg) : m_cfg(cfg) {
  Stepper stepper(m_cfg.magneticField);
  Navigator::Config navCfg;
  navCfg.detector = m_cfg.detector;
  navCfg.resolvePassive = false;
  navCfg.resolveMaterial = true;
  navCfg.resolveSensitive = true;
  Navigator navigator(navCfg,
                      Acts::getDefaultLogger("AlignmentDetectorNavigator",
                                             Acts::Logging::INFO));
  Propagator propagator(std::move(stepper), std::move(navigator));
  KFFitter trackFitter(
      std::move(propagator),
      Acts::getDefaultLogger("AlignmentKalmanFilter", Acts::Logging::INFO));
  m_align = std::make_unique<ActsAligner>(std::move(trackFitter));
}

AlignmentFunction::AlignmentResult ActsAlignmentFunction::operator()(
    const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
    const Acts::CalibrationContext& cctx,
    const SourceLinkContainer& alignmentFitSourceLinks,
    const SourceLinkContainer& initialTrackStateFitSourceLinks,
    const TrackParametersContainer& initialParameters,
    const MagneticFieldParametersContainer& magFieldParameters) const {
  // Initialize KF options
  auto propOptions = KFFitterPropagatorOptions(gctx, mctx);
  propOptions.maxSteps = m_cfg.maxKFSteps;

  auto kfOptions = Acts::KalmanFitterOptions(gctx, mctx, cctx,
                                             m_cfg.kfExtensions, propOptions);

  kfOptions.referenceSurface = m_cfg.kfReferenceSurface;

  // Set the alignment options
  ActsAlignment::AlignmentOptions<KFFitterOptions> alignOptions(
      kfOptions, m_cfg.alignmentTransformUpdater,
      m_cfg.alignmentParametersSolver, m_cfg.alignedDetElements,
      m_cfg.chi2ONdfCutOff, m_cfg.deltaChi2ONdfCutOff,
      m_cfg.maxAlignmentFitNumIt, m_cfg.alignmentMask);

  // Alignment fit loop
  std::vector<Acts::GenericCurvilinearTrackParameters<Acts::ParticleHypothesis>>
      trackParameters(initialParameters.begin(), initialParameters.end());
  ActsAlignmentResult actsAlignmentResult;
  for (std::size_t i = 0; i < m_cfg.nRefittingIt; i++) {
    // Run alignment
    actsAlignmentResult = m_align
                              ->align(alignmentFitSourceLinks, trackParameters,
                                      magFieldParameters, alignOptions)
                              .value();

    // Re-estimate the track parameters
    for (std::size_t j = 0; j < initialTrackStateFitSourceLinks.size(); j++) {
      const auto& candidateSourceLinks = initialTrackStateFitSourceLinks.at(j);
      const auto& originParameters = trackParameters.at(j);
      Acts::MagneticFieldContext mctxCurrent(magFieldParameters.at(j));

      trackParameters.insert(
          trackParameters.begin() + j,
          m_cfg.trackParametersEstimator
              ->estimateParameters(gctx, mctxCurrent, candidateSourceLinks,
                                   originParameters.direction(),
                                   originParameters.position(gctx))
              .value());
    }
  }

  // Fill out alignment results struct
  AlignmentResult alignmentResult;
  alignmentResult.deltaAlignmentParameters =
      actsAlignmentResult.deltaAlignmentParameters;
  alignmentResult.alignedParameters = actsAlignmentResult.alignedParameters;
  alignmentResult.alignmentCovariance = actsAlignmentResult.alignmentCovariance;
  alignmentResult.averageChi2ONdf = actsAlignmentResult.averageChi2ONdf;
  alignmentResult.deltaChi2 = actsAlignmentResult.deltaChi2;
  alignmentResult.chi2 = actsAlignmentResult.chi2;
  alignmentResult.measurementDim = actsAlignmentResult.measurementDim;
  alignmentResult.alignmentDof = actsAlignmentResult.alignmentDof;
  alignmentResult.numTracks = actsAlignmentResult.numTracks;
  alignmentResult.updatedInitialTrackStates = trackParameters;

  return alignmentResult;
}
