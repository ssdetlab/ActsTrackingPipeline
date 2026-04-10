#pragma once

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/AbortList.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"

/// -----------------------------------------------
/// Propagator definitions

using Navigator = Acts::Experimental::DetectorNavigator;
using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Navigator>;

/// -----------------------------------------------
/// KF fitter definitions

using KFFitterGainUpdater = Acts::GainMatrixUpdater;
using KFFitterGainSmoother = Acts::GainMatrixSmoother;

using KFFitterActionList = Acts::ActionList<>;
using KFFitterAbortList = Acts::AbortList<Acts::EndOfWorldReached>;

using KFFitterPropagator =
    Acts::Propagator<Acts::EigenStepper<>,
                     Acts::Experimental::DetectorNavigator>;

using KFFitterTrackContainerBackend = Acts::VectorTrackContainer;
using KFFitterTrajectory = Acts::VectorMultiTrajectory;

using KFFitterOptions = Acts::KalmanFitterOptions<KFFitterTrajectory>;

using KFFitterTrackContainer =
    Acts::TrackContainer<KFFitterTrackContainerBackend, KFFitterTrajectory,
                         std::shared_ptr>;

using KFFitterPropagatorOptions =
    typename Propagator::template Options<KFFitterActionList,
                                          KFFitterAbortList>;

using KFFitterExtensions = Acts::KalmanFitterExtensions<KFFitterTrajectory>;

using KFFitter = Acts::KalmanFitter<Propagator, KFFitterTrajectory>;

/// -----------------------------------------------
/// Alignment fitter definitions

using Alignment = ActsAlignment::Alignment<KFFitter>;
using AlignmentResult = ActsAlignment::AlignmentResult;
