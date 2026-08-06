#include "TrackingPipeline/Alignment/MinuitGlobalAlignmentFunction.hpp"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/detail/AlignmentStoreBuilders.hpp"
#include "TrackingPipeline/EventData/ExtendedSourceLink.hpp"

MinuitGlobalAlignmentFunction::MinuitAlignmentFunctionImpl::
    MinuitAlignmentFunctionImpl(const Options& opt)
    : m_opt(opt) {}

double MinuitGlobalAlignmentFunction::MinuitAlignmentFunctionImpl::Up() const {
  return 1;
}

double MinuitGlobalAlignmentFunction::MinuitAlignmentFunctionImpl::operator()(
    const std::vector<double>& pars) const {
  Acts::GeometryContext gctxDef;
  Acts::MagneticFieldContext mctxDef;
  Acts::CalibrationContext cctxDef;

  std::cout << "--------------------------------------------------\n";
  std::cout << "OPERATOR ()\n";
  std::cout << "PARS:\n";
  for (double par : pars) {
    std::cout << par << "\n";
  }

  double dy = pars.at(0);
  double dz = pars.at(1);
  double dtheta = pars.at(2);
  std::cout << "DY " << dy << "\n";
  std::cout << "DZ " << dz << "\n";
  std::cout << "DTHETA " << dtheta << "\n";

  // Initialize alignment store
  Acts::Vector3 globalShiftMean(0, dy, dz);
  std::unordered_map<int, Acts::Vector3> localShiftsMean{
      {10, Acts::Vector3(0, 0, 0)},
      {12, Acts::Vector3(0, 0, 0)},
      {14, Acts::Vector3(0, 0, 0)},
      {16, Acts::Vector3(0, 0, 0)},
      {18, Acts::Vector3(0, 0, 0)}};
  Acts::Vector3 globalAnglesMean(0, 0, dtheta);
  std::unordered_map<int, Acts::Vector3> localAnglesMean{
      {10, Acts::Vector3(0, 0, 0)},
      {12, Acts::Vector3(0, 0, 0)},
      {14, Acts::Vector3(0, 0, 0)},
      {16, Acts::Vector3(0, 0, 0)},
      {18, Acts::Vector3(0, 0, 0)}};

  auto aStore = detail::makeAlignmentStore(gctxDef, m_opt.detector,
                                           globalShiftMean, localShiftsMean,
                                           globalAnglesMean, localAnglesMean);
  AlignmentContext alignCtx(aStore);
  Acts::GeometryContext gctx = Acts::GeometryContext{alignCtx};

  // Estimate the track parameters
  std::vector<Acts::GenericCurvilinearTrackParameters<Acts::ParticleHypothesis>>
      trackParameters;
  trackParameters.reserve(m_opt.trackParametersFitSourceLinkContainer.size());
  for (std::size_t j = 0;
       j < m_opt.trackParametersFitSourceLinkContainer.size(); j++) {
    const auto& candidateSourceLinks =
        m_opt.trackParametersFitSourceLinkContainer.at(j);
    const auto& originParameters = m_opt.trackParameters.at(j);
    Acts::MagneticFieldContext mctxCurrent(m_opt.magFieldParameters.at(j));

    trackParameters.push_back(
        m_opt.trackParametersEstimator
            ->estimateParameters(gctx, mctxCurrent, candidateSourceLinks,
                                 originParameters.direction(),
                                 originParameters.position(gctx))
            .value());
  }
  std::cout << "NEW TRACK PARS " << trackParameters.size() << "\n";

  // Initialize KF options
  auto propOptions = KFFitterPropagatorOptions(gctx, mctxDef);
  propOptions.maxSteps = m_opt.maxKFSteps;

  auto kfOptions =
      Acts::KalmanFitterOptions(gctx, mctxDef, cctxDef, m_opt.kfExtensions,
                                propOptions, m_opt.kfReferenceSurface);

  double chi2Smoothed = 0;
  std::cout << "STARTING FIT\n";
  for (unsigned int iTraj = 0; iTraj < m_opt.trackFitSourceLinkContainer.size();
       iTraj++) {
    const auto& sls = m_opt.trackFitSourceLinkContainer.at(iTraj);
    const auto& sParameters = trackParameters.at(iTraj);
    const auto& mParameters = m_opt.magFieldParameters.at(iTraj);

    Acts::MagneticFieldContext trajMctx{mParameters};
    kfOptions.magFieldContext =
        std::reference_wrapper<const Acts::MagneticFieldContext>(trajMctx);

    Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                                Acts::VectorMultiTrajectory{}};
    // Perform the fit
    auto fitRes = m_opt.fitter.fit(sls.begin(), sls.end(), sParameters,
                                   kfOptions, tracks);
    if (!fitRes.ok()) {
      std::cout << "Fit failure\n";
      continue;
    }
    // The fit results
    const auto& track = fitRes.value();
    for (const auto& state : track.trackStatesReversed()) {
      // Skip the states without meaningful information
      if (!state.hasProjector()) {
        continue;
      }
      Acts::ActsDynamicVector meas = state.effectiveCalibrated();
      Acts::ActsDynamicMatrix measurementCov =
          state.effectiveCalibratedCovariance();

      Acts::ActsDynamicMatrix smoothedState =
          state.effectiveProjector() * state.smoothed();
      Acts::ActsDynamicMatrix smoothedCov =
          measurementCov - state.effectiveProjector() *
                               state.smoothedCovariance() *
                               state.effectiveProjector().transpose();
      Acts::ActsDynamicVector smoothedDiag =
          smoothedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();

      Acts::ActsDynamicVector smoothedPull =
          smoothedDiag.cwiseProduct(meas - smoothedState);

      chi2Smoothed += smoothedPull.dot(smoothedPull);
    }
    if (iTraj % 5000 == 0) {
      std::cout << "TRACK E: " << std::abs(1.0 / track.qOverP()) << "\n";
    }
  }

  std::cout << "CHI2 PER TRACK "
            << chi2Smoothed / m_opt.trackFitSourceLinkContainer.size() << "\n";
  std::cout << "--------------------------------------------------\n";
  int dummy = 0;
  return chi2Smoothed;
}

MinuitGlobalAlignmentFunction::MinuitGlobalAlignmentFunction(const Config& cfg)
    : m_cfg(cfg) {}

AlignmentResult MinuitGlobalAlignmentFunction::operator()(
    const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
    const Acts::CalibrationContext& cctx,
    const AlignmentAlgorithm::SourceLinkContainer& sourceLinks,
    const AlignmentAlgorithm::TrackParametersContainer& initialParameters,
    const AlignmentAlgorithm::MagneticFieldParametersContainer&
        magFieldParameters) {
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

  auto typeRange = sourceLinks.front().getTypeRange();
  std::vector<std::vector<std::size_t>> idxRanges;

  AlignmentAlgorithm::SourceLinkContainer tfSourceLinks;
  tfSourceLinks.reserve(sourceLinks.size());
  for (const auto& slc : sourceLinks) {
    // std::cout << "------------------------------------------------------------\n";
    std::vector<std::size_t> sourceLinksIndices = slc.getIndexRange();
    std::vector<std::size_t> filteredSourceLinksIndices;
    filteredSourceLinksIndices.reserve(5);
    // std::cout << "INITIAL: " << sourceLinksIndices.size() << "\n";
    for (std::size_t i = 0; i < sourceLinksIndices.size(); i++) {
      std::size_t idx = sourceLinksIndices.at(i);
      // std::cout << "IDX: " << i << " -- " << idx << "\n";
      if (slc.at(i).type() == typeid(ExtendedSourceLink)) {
        continue;
      }
      // std::cout << "ACCEPTED\n";
      filteredSourceLinksIndices.push_back(idx);
    }
    // std::cout << "FILTERED: " << filteredSourceLinksIndices.size() << "\n";
    // --------------------------------------------------

    idxRanges.push_back(filteredSourceLinksIndices);
  }
  for (const auto& range : idxRanges) {
    tfSourceLinks.emplace_back(typeRange, range);
  }

  MinuitAlignmentFunctionImpl::Options alignmentFunctionCfg{
      .trackFitSourceLinkContainer = sourceLinks,
      .trackParametersFitSourceLinkContainer = tfSourceLinks,
      .trackParameters = initialParameters,
      .magFieldParameters = magFieldParameters,
      .fitter = trackFitter,
      .kfExtensions = m_cfg.kfExtensions,
      .kfReferenceSurface = m_cfg.kfReferenceSurface,
      .maxKFSteps = m_cfg.maxKFSteps,
      .trackParametersEstimator = m_cfg.trackParametersEstimator,
      .detector = m_cfg.detector};
  MinuitAlignmentFunctionImpl alignmentFunction(alignmentFunctionCfg);

  // Create Minuit parameters with names
  Acts::Vector3 pars(10, -5, -2e-3);
  Acts::Vector3 parsErrs = Acts::Vector3(1, 1, 1e-3);
  ROOT::Minuit2::MnUserParameters userPars;
  userPars.Add("deltaGlobalY", pars(0), parsErrs(0), -20, 20);
  userPars.Add("deltaGlobalZ", pars(1), parsErrs(1), -20, 20);
  userPars.Add("deltaGlobalTheta", pars(2), parsErrs(2), -5e-3, 5e-3);

  // Create MIGRAD minimizer
  ROOT::Minuit2::MnMigrad migrad(alignmentFunction, userPars);
  // Minimize
  ROOT::Minuit2::FunctionMinimum min = migrad();
  throw std::runtime_error("ERR");
}
