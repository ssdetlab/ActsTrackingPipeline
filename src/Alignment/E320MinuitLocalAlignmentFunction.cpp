#include "TrackingPipeline/Alignment/E320MinuitLocalAlignmentFunction.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <memory>

#include <Eigen/Core>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnHesse.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnUserParameters.h"
#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/AlignmentFunction.hpp"

E320::E320MinuitLocalAlignmentFunction::MinuitAlignmentFunctionImpl::
    MinuitAlignmentFunctionImpl(const Options& opt, Acts::Logging::Level level)
    : m_opt(opt),
      m_logger(Acts::getDefaultLogger("MinuitAlignmentFunctionImpl", level)) {
  // Initialize internal state
  m_state = std::make_unique<State>();

  // Compute detector COM
  m_com = Acts::Vector3::Zero();
  for (const auto* element : m_opt.functionCfg.alignedDetElements) {
    const auto& surf = element->surface();
    m_com += surf.center(m_opt.gctx);
  }
  m_com /= m_opt.functionCfg.alignedDetElements.size();

  for (const auto& element : m_opt.functionCfg.alignedDetElements) {
    const auto& surf = element->surface();
    m_leverArms[surf.geometryId()] = surf.center(m_opt.gctx) - m_com;
  }

  // Store nominal alignment parameters
  auto nominalAlignmentStore = *m_opt.functionCfg.alignmentStore;
  m_nominalGctx = Acts::GeometryContext(AlignmentContext(
      std::make_shared<AlignmentStore>(nominalAlignmentStore)));
}

const E320::E320MinuitLocalAlignmentFunction::MinuitAlignmentFunctionImpl::
    State&
    E320::E320MinuitLocalAlignmentFunction::MinuitAlignmentFunctionImpl::state()
        const {
  return *m_state;
}

double E320::E320MinuitLocalAlignmentFunction::MinuitAlignmentFunctionImpl::Up()
    const {
  return m_opt.functionCfg.upLevel;
}

double
E320::E320MinuitLocalAlignmentFunction::MinuitAlignmentFunctionImpl::operator()(
    const std::vector<double>& pars) const {
  // Initialize KF options
  auto propOptions = KFFitterPropagatorOptions(m_opt.gctx, m_opt.mctx);
  propOptions.maxSteps = m_opt.functionCfg.maxKFSteps;

  auto kfOptions = Acts::KalmanFitterOptions(
      m_opt.gctx, m_opt.mctx, m_opt.cctx, m_opt.functionCfg.kfExtensions,
      propOptions, m_opt.functionCfg.kfReferenceSurface);

  ACTS_INFO("-----------------------------------------");
  Acts::ActsDynamicVector deltas(pars.size());
  for (std::size_t i = 0; i < pars.size(); i++) {
    deltas(i) = pars.at(i);
    ACTS_INFO("Parameter " << i << ": " << pars.at(i));
  }
  std::cout << "DELTAS: " << deltas.transpose() << "\n";

  // Update alignment parameters
  for (std::size_t i = 0; i < m_opt.functionCfg.alignedDetElements.size();
       i++) {
    const auto* element = m_opt.functionCfg.alignedDetElements.at(i);
    const auto& surf = element->surface();
    const Acts::Transform3& nominalTransform = surf.transform(m_nominalGctx);
    const Acts::RotationMatrix3& nominalRotation = nominalTransform.rotation();

    Acts::Vector3 deltaPars = deltas.block(i * 3, 0, 3, 1);
    std::cout << "DELTA PARS " << deltaPars.transpose() << "\n";

    Acts::RotationMatrix3 deltaRotation =
        Acts::AngleAxis3(deltaPars(2), Acts::Vector3::UnitZ())
            .toRotationMatrix() *
        Acts::AngleAxis3(0, Acts::Vector3::UnitY()).toRotationMatrix() *
        Acts::AngleAxis3(0, Acts::Vector3::UnitX()).toRotationMatrix();

    Acts::Vector3 deltaTranslation =
        Acts::Vector3(0, deltaPars(0), deltaPars(1));

    Acts::Transform3 newTransform = Acts::Transform3::Identity();
    Acts::Vector3 newCenter =
        m_com + m_leverArms.at(surf.geometryId()) + deltaTranslation;
    Acts::RotationMatrix3 newRotation = nominalRotation * deltaRotation;

    newTransform.translation() = newCenter;
    newTransform.rotate(newRotation);
    std::cout << "NEW CENTER " << newTransform.translation().transpose()
              << "\n";

    m_opt.functionCfg.alignmentStore->store[element->surface().geometryId()] =
        newTransform;

    ACTS_INFO("-----------------------------------------");
    ACTS_INFO("Surface " << surf.geometryId());
    ACTS_INFO("Delta angles\n" << deltaRotation);
    ACTS_INFO("Delta translation " << deltaTranslation.transpose());
    ACTS_INFO("Center " << surf.center(m_opt.gctx).transpose());
    ACTS_INFO("Normal " << surf.normal(m_opt.gctx, surf.center(m_opt.gctx),
                                       Acts::Vector3::UnitX())
                               .transpose());
    ACTS_INFO("Rotation \n" << surf.transform(m_opt.gctx).rotation());
  }

  // Estimate the track parameters
  std::vector<Acts::GenericCurvilinearTrackParameters<Acts::ParticleHypothesis>>
      trackParameters;
  trackParameters.reserve(m_opt.initialTrackStateFitSourceLinkContainer.size());
  for (std::size_t j = 0;
       j < m_opt.initialTrackStateFitSourceLinkContainer.size(); j++) {
    const auto& candidateSourceLinks =
        m_opt.initialTrackStateFitSourceLinkContainer.at(j);
    const auto& originParameters = m_opt.trackParameters.at(j);
    Acts::MagneticFieldContext mctxCurrent(m_opt.magFieldParameters.at(j));

    trackParameters.push_back(
        m_opt.functionCfg.trackParametersEstimator
            ->estimateParameters(m_opt.gctx, mctxCurrent, candidateSourceLinks,
                                 originParameters.direction(),
                                 originParameters.position(m_opt.gctx))
            .value());
  }

  // Run the KF fit
  std::size_t nTracks = m_opt.alignmentFitSourceLinkContainer.size();

  double chi2Smoothed = 0;
  std::size_t measurementDim = 0;

  Acts::SquareMatrix4 matEXY = Acts::SquareMatrix4::Zero();
  Acts::Vector4 vecEX = Acts::Vector4::Zero();

  std::vector<Acts::Vector4> ipRes;
  ipRes.reserve(nTracks);

  Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                              Acts::VectorMultiTrajectory{}};
  tracks.container().reserve(nTracks);
  for (std::size_t iTraj = 0; iTraj < nTracks; iTraj++) {
    const auto& sls = m_opt.alignmentFitSourceLinkContainer.at(iTraj);
    const auto& sParameters = trackParameters.at(iTraj);
    const auto& mParameters = m_opt.magFieldParameters.at(iTraj);

    Acts::MagneticFieldContext trajMctx{mParameters};
    kfOptions.magFieldContext =
        std::reference_wrapper<const Acts::MagneticFieldContext>(trajMctx);

    // Perform the track fit
    auto fitRes = m_opt.fitter->fit(sls.begin(), sls.end(), sParameters,
                                    kfOptions, tracks);
    if (!fitRes.ok()) {
      std::cout << "Fit failure\n";
      continue;
    }

    // Fit results
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

      measurementDim += state.calibratedSize();
    }
  }
  ACTS_INFO("Chi2 per track:" << chi2Smoothed / nTracks);

  // Update the function state
  m_state->chi2 = chi2Smoothed;
  m_state->measurementDim = measurementDim;

  return chi2Smoothed;
}

E320::E320MinuitLocalAlignmentFunction::E320MinuitLocalAlignmentFunction(
    const Config& cfg, Acts::Logging::Level level)
    : m_cfg(cfg),
      m_level(level),
      m_logger(
          Acts::getDefaultLogger("E320MinuitLocalAlignmentFunction", level)) {
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
  m_fitter = std::make_shared<KFFitter>(
      std::move(propagator),
      Acts::getDefaultLogger("AlignmentKalmanFilter", Acts::Logging::INFO));
}

AlignmentFunction::AlignmentResult
E320::E320MinuitLocalAlignmentFunction::operator()(
    const Acts::GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
    const Acts::CalibrationContext& cctx,
    const SourceLinkContainer& alignmentFitSourceLinks,
    const SourceLinkContainer& initialTrackStateFitSourceLinks,
    const TrackParametersContainer& initialParameters,
    const MagneticFieldParametersContainer& magFieldParameters) const {
  // Store nominal alignment parameters
  auto nominalAlignmentStore = *m_cfg.alignmentStore;
  Acts::GeometryContext nominalGctx(AlignmentContext(
      std::make_shared<AlignmentStore>(nominalAlignmentStore)));

  // Set alignment options
  MinuitAlignmentFunctionImpl::Options alignmentFunctionCfg{
      .gctx = gctx,
      .mctx = mctx,
      .cctx = cctx,
      .alignmentFitSourceLinkContainer = alignmentFitSourceLinks,
      .initialTrackStateFitSourceLinkContainer =
          initialTrackStateFitSourceLinks,
      .trackParameters = initialParameters,
      .magFieldParameters = magFieldParameters,
      .fitter = m_fitter,
      .functionCfg = m_cfg};
  MinuitAlignmentFunctionImpl alignmentFunction(alignmentFunctionCfg, m_level);

  // Evaluate the FCN to get the initial state
  alignmentFunction(
      {m_cfg.initialPars(0), m_cfg.initialPars(1), m_cfg.initialPars(2),
       m_cfg.initialPars(3), m_cfg.initialPars(4), m_cfg.initialPars(5),
       m_cfg.initialPars(6), m_cfg.initialPars(7), m_cfg.initialPars(8),
       m_cfg.initialPars(9), m_cfg.initialPars(10), m_cfg.initialPars(11)});
  const auto& state = alignmentFunction.state();
  double initialChi2 = state.chi2;
  std::size_t measurementDim = state.measurementDim;

  // Create Minuit parameters with names
  ROOT::Minuit2::MnUserParameters userPars;
  userPars.Add("dy0", m_cfg.initialPars(0), m_cfg.initialParSteps(0),
               m_cfg.parLowBounds(0), m_cfg.parHighBounds(0));
  userPars.Add("dz0", m_cfg.initialPars(1), m_cfg.initialParSteps(1),
               m_cfg.parLowBounds(1), m_cfg.parHighBounds(1));
  userPars.Add("dTheta0", m_cfg.initialPars(2), m_cfg.initialParSteps(2),
               m_cfg.parLowBounds(2), m_cfg.parHighBounds(2));

  userPars.Add("dy1", m_cfg.initialPars(3), m_cfg.initialParSteps(3),
               m_cfg.parLowBounds(3), m_cfg.parHighBounds(3));
  userPars.Add("dz1", m_cfg.initialPars(4), m_cfg.initialParSteps(4),
               m_cfg.parLowBounds(4), m_cfg.parHighBounds(4));
  userPars.Add("dTheta1", m_cfg.initialPars(5), m_cfg.initialParSteps(5),
               m_cfg.parLowBounds(5), m_cfg.parHighBounds(5));

  userPars.Add("dy2", m_cfg.initialPars(6), m_cfg.initialParSteps(6),
               m_cfg.parLowBounds(6), m_cfg.parHighBounds(6));
  userPars.Add("dz2", m_cfg.initialPars(7), m_cfg.initialParSteps(7),
               m_cfg.parLowBounds(7), m_cfg.parHighBounds(7));
  userPars.Add("dTheta2", m_cfg.initialPars(8), m_cfg.initialParSteps(8),
               m_cfg.parLowBounds(8), m_cfg.parHighBounds(8));

  userPars.Add("dy3", m_cfg.initialPars(9), m_cfg.initialParSteps(9),
               m_cfg.parLowBounds(9), m_cfg.parHighBounds(9));
  userPars.Add("dz3", m_cfg.initialPars(10), m_cfg.initialParSteps(10),
               m_cfg.parLowBounds(10), m_cfg.parHighBounds(10));
  userPars.Add("dTheta3", m_cfg.initialPars(11), m_cfg.initialParSteps(11),
               m_cfg.parLowBounds(11), m_cfg.parHighBounds(11));

  // Create MIGRAD minimizer
  ROOT::Minuit2::MnMigrad migrad(alignmentFunction, userPars);
  // Minimize
  ROOT::Minuit2::FunctionMinimum min = migrad();

  // Create HESSE estimator
  ROOT::Minuit2::MnHesse hesse;
  hesse(alignmentFunction, min);

  // Output
  ACTS_INFO("Finished MIGRAD minimization with chi2: " << min.Fval());
  ACTS_INFO("HESSE min: " << min);

  // Estimate the track parameters
  std::vector<Acts::GenericCurvilinearTrackParameters<Acts::ParticleHypothesis>>
      trackParameters;
  trackParameters.reserve(initialTrackStateFitSourceLinks.size());
  for (std::size_t j = 0; j < initialTrackStateFitSourceLinks.size(); j++) {
    const auto& candidateSourceLinks = initialTrackStateFitSourceLinks.at(j);
    const auto& originParameters = initialParameters.at(j);
    Acts::MagneticFieldContext mctxCurrent(magFieldParameters.at(j));

    trackParameters.push_back(
        m_cfg.trackParametersEstimator
            ->estimateParameters(gctx, mctxCurrent, candidateSourceLinks,
                                 originParameters.direction(),
                                 originParameters.position(gctx))
            .value());
  }

  // Store aligned transforms
  std::unordered_map<Acts::DetectorElementBase*, Acts::Transform3>
      alignedParameters;
  for (auto* element : m_cfg.alignedDetElements) {
    alignedParameters.emplace(element, element->surface().transform(gctx));
  }

  // Fill out alignment results struct
  double finalChi2 = min.Fval();
  std::size_t ndf = 0;

  // TODO: Finish the alignment result
  AlignmentResult alignmentResult;
  alignmentResult.deltaAlignmentParameters = Acts::Vector3(
      min.UserParameters().Value(0), min.UserParameters().Value(1),
      min.UserParameters().Value(2));
  alignmentResult.alignedParameters = alignedParameters;
  alignmentResult.averageChi2ONdf = finalChi2 / measurementDim;
  alignmentResult.deltaChi2 = initialChi2 - finalChi2;
  alignmentResult.chi2 = finalChi2;
  alignmentResult.measurementDim = measurementDim;
  alignmentResult.alignmentDof = userPars.Parameters().size();
  alignmentResult.numTracks = initialParameters.size();
  alignmentResult.updatedInitialTrackStates = trackParameters;

  std::size_t nRows = min.UserCovariance().Nrow();
  alignmentResult.alignmentCovariance = Acts::ActsDynamicMatrix(nRows, nRows);
  for (std::size_t i = 0; i < nRows; i++) {
    for (std::size_t j = 0; j < nRows; j++) {
      alignmentResult.alignmentCovariance(i, j) = min.UserCovariance()(i, j);
    }
  }

  return alignmentResult;
}
