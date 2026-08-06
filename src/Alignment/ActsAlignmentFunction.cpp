#include "TrackingPipeline/Alignment/ActsAlignmentFunction.hpp"

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
  m_align = std::make_unique<Alignment>(std::move(trackFitter));
}

AlignmentResult ActsAlignmentFunction::operator()(
    const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
    const Acts::CalibrationContext& cctx,
    const AlignmentAlgorithm::SourceLinkContainer& sourceLinks,
    const AlignmentAlgorithm::TrackParametersContainer& initialParameters,
    const AlignmentAlgorithm::MagneticFieldParametersContainer&
        magFieldParameters) {
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

  // Run alignment
  return m_align
      ->align(sourceLinks, initialParameters, magFieldParameters, alignOptions)
      .value();
}
