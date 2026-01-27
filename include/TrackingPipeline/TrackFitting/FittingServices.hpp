#pragma once

#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

struct FittingServices {
  using Stepper    = Acts::EigenStepper<>;
  using Navigator  = Acts::Experimental::DetectorNavigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using Trajectory = Acts::VectorMultiTrajectory;

  std::shared_ptr<const Acts::Experimental::Detector> detector;

  std::shared_ptr<Propagator> propagator;
  std::shared_ptr<Acts::KalmanFitter<Propagator, Trajectory>> kalmanFitter;
  std::optional<Acts::KalmanFitterOptions<Trajectory>> baseKfOptions;
  std::optional<Acts::KalmanFitterExtensions<Trajectory>> baseExtensions;

  // for writers
  std::optional<SimpleSourceLink::SurfaceAccessor> surfaceAccessor;
  std::shared_ptr<Acts::Surface> referenceSurface;

  static FittingServices& instance() {
    static FittingServices s;
    return s;
  }
};
