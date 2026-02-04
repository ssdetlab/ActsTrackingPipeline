#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Simulation/SimpleDigitizer.hpp"
#include "TrackingPipeline/Simulation/GaussianVertexGenerator.hpp"
#include "TrackingPipeline/Simulation/SphericalMomentumGenerator.hpp"

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

  // --- FastSim helpers ---
  std::shared_ptr<SimpleDigitizer>              simDigitizer;
  std::shared_ptr<GaussianVertexGenerator>      simVertexGenerator;
  std::shared_ptr<SphericalMomentumGenerator>   simMomentumGenerator;
  // Sensitive detector surfaces for background generator
  std::vector<const Acts::Surface*> simDetSurfaces;

  // --- Lookup helpers ---
  // Reference tracking layers for lookup-table estimation
  std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*>
      lookupRefLayers;

  static FittingServices& instance() {
    static FittingServices s;
    return s;
  }
};
