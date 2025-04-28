#pragma once

#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <ActsAlignment/Kernel/Alignment.hpp>

#include <cstddef>
#include <memory>

#include "TMinuit.h"
#include "TrackingPipeline/Alignment/AlignmentAlgorithm.hpp"

namespace MinuitAlignment {

using Updater = Acts::GainMatrixUpdater;
using Smoother = Acts::GainMatrixSmoother;
using Stepper = Acts::EigenStepper<>;
using Propagator =
    Acts::Propagator<Stepper, Acts::Experimental::DetectorNavigator>;
using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;

class MinuitAlignment {
 public:
  MinuitAlignment(
      std::shared_ptr<const Acts::Experimental::Detector> detector,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField) {
    Stepper stepper(std::move(magneticField));
    Acts::Experimental::DetectorNavigator::Config cfg;
    cfg.detector = detector.get();
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    Acts::Experimental::DetectorNavigator navigator(
        cfg, Acts::getDefaultLogger("DetectorNavigator", Acts::Logging::INFO));
    Propagator propagator(std::move(stepper), std::move(navigator));
    m_fitter = std::make_shared<Fitter>(std::move(propagator));

    m_minuit = new TMinuit();
    m_minuit->SetFCN(costFunc);

    m_minuit->mnparm(0, "dx0", 0, 1e-2, -10, 10, m_ierflag);
    /*m_minuit->mnparm(1, "dx1", 0, 1e-2, -10, 10, m_ierflag);*/
    /*m_minuit->mnparm(2, "dy0", 0, 1e-2, -10, 10, m_ierflag);*/
    /*m_minuit->mnparm(3, "dy1", 0, 1e-2, -10, 10, m_ierflag);*/

    m_trajectories = std::vector<std::vector<Acts::SourceLink>>{};
    m_initialParameters = std::vector<Acts::CurvilinearTrackParameters>{};
    m_options = nullptr;
  }

  Acts::Result<ActsAlignment::AlignmentResult> align(
      const std::vector<std::vector<Acts::SourceLink>>& sourceLinks,
      const std::vector<Acts::CurvilinearTrackParameters>& initialParameters,
      const ActsAlignment::AlignmentOptions<
          AlignmentAlgorithm::TrackFitterOptions>& options) {
    for (unsigned int iDetElement = 0;
         iDetElement < options.alignedDetElements.size(); iDetElement++) {
      m_result.idxedAlignSurfaces.emplace(
          &options.alignedDetElements.at(iDetElement)->surface(), iDetElement);
    }
    std::cout << "There are " << m_result.idxedAlignSurfaces.size()
              << " detector elements to be aligned\n";

    std::cout << "Max number of iterations: " << options.maxIterations << "\n";
    double chi2 = 0;
    m_trajectories = sourceLinks;
    m_initialParameters = initialParameters;
    m_options = std::make_shared<ActsAlignment::AlignmentOptions<
        AlignmentAlgorithm::TrackFitterOptions>>(options);

    m_oldTransform0 = m_options->alignedDetElements.at(0)->transform(
        m_options->fitOptions.geoContext);
    /*m_oldTransform1 = m_options->alignedDetElements.at(1)->transform(*/
    /*    m_options->fitOptions.geoContext);*/

    m_minuit->Command("SET ERR 1");
    m_minuit->Command("SET PRINT 2");
    m_minuit->Command("SET STR 1");
    m_minuit->Command("MIGRAD 1000 0.01");
    /*m_minuit->Command("IMPROVE 1000");*/

    return m_result;
  }

  static void costFunc(int& npar, double* gin, double& f, double* par,
                       int iflag) {
    double chi2 = 0;

    m_idx++;
    double x0 = par[0];
    /*double x1 = par[1];*/
    /*double y0 = par[2];*/
    /*double y1 = par[3];*/
    std::cout << "dx0 = " << x0 << "\n";
    /*std::cout << "dx1 = " << x1 << "\n";*/
    /*std::cout << "dy0 = " << y0 << "\n";*/
    /*std::cout << "dy1 = " << y1 << "\n";*/

    Acts::Vector3 newCenter0 =
        Acts::Vector3(m_oldTransform0.translation().x() - x0,
                      m_oldTransform0.translation().y(),
                      m_oldTransform0.translation().z());
    /*Acts::Vector3 newCenter1 =*/
    /*    Acts::Vector3(m_oldTransform1.translation().x() - x1,*/
    /*                  m_oldTransform1.translation().y(),*/
    /*                  m_oldTransform1.translation().z() - y1);*/

    Acts::Transform3 newTransform0 = m_oldTransform0;
    newTransform0.translation() = newCenter0;

    /*Acts::Transform3 newTransform1 = m_oldTransform1;*/
    /*newTransform1.translation() = newCenter1;*/

    bool updated1 = m_options->alignedTransformUpdater(
        m_options->alignedDetElements.at(0), m_options->fitOptions.geoContext,
        newTransform0);
    /*bool updated0 = m_options->alignedTransformUpdater(*/
    /*    m_options->alignedDetElements.at(1), m_options->fitOptions.geoContext,*/
    /*    newTransform1);*/

    for (std::size_t iTraj = 0; iTraj < m_trajectories.size(); iTraj++) {
      const auto& sourcelinks = m_trajectories.at(iTraj);
      const auto& sParameters = m_initialParameters.at(iTraj);
      // The result for one single track
      Acts::TrackContainer tracks{Acts::VectorTrackContainer{},
                                  Acts::VectorMultiTrajectory{}};

      // Convert to Acts::SourceLink during iteration
      Acts::SourceLinkAdapterIterator begin{sourcelinks.begin()};
      Acts::SourceLinkAdapterIterator end{sourcelinks.end()};

      // Perform the fit
      auto fitRes =
          m_fitter->fit(begin, end, sParameters, m_options->fitOptions, tracks);

      if (!fitRes.ok()) {
        std::cout << "Fit failure\n";
      }
      // The fit results
      const auto& track = fitRes.value();
      for (const auto& state : track.trackStatesReversed()) {
        if (!state.hasProjector()) {
          continue;
        }
        // Get the measurements hit
        auto hit = state.effectiveCalibrated();
        auto smoothedHit = state.effectiveProjector() * state.smoothed();
        auto measurementCov = state.effectiveCalibratedCovariance();
        auto smoothedCov = state.effectiveProjector() *
                               state.smoothedCovariance() *
                               state.effectiveProjector().transpose() -
                           measurementCov;
        auto smoothedDiag =
            smoothedCov.cwiseAbs().diagonal().cwiseInverse().cwiseSqrt();
        auto smoothedPull = smoothedDiag.cwiseProduct(hit - smoothedHit);
        chi2 += smoothedPull.dot(smoothedPull);
      }
    }
    std::cout << "CHI2 " << chi2 << "\n";
    f = chi2;
  }

 private:
  // The fitter
  inline static std::shared_ptr<Fitter> m_fitter;

  inline static int m_ierflag = 0;
  inline static std::size_t m_idx = 0;
  TMinuit* m_minuit = nullptr;

  inline static Acts::Transform3 m_oldTransform0;
  inline static Acts::Transform3 m_oldTransform1;

  inline static std::vector<std::vector<Acts::SourceLink>> m_trajectories;
  inline static std::vector<Acts::CurvilinearTrackParameters>
      m_initialParameters;
  inline static std::shared_ptr<
      ActsAlignment::AlignmentOptions<AlignmentAlgorithm::TrackFitterOptions>>
      m_options;
  ActsAlignment::AlignmentResult m_result;
};

}  // namespace MinuitAlignment
