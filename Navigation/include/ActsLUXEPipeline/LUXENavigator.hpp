#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/GenericCurvilinearTrackParameters.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"


#include <iostream>
#include <vector>
#include <functional>

namespace LUXENavigator {
using namespace Acts::UnitLiterals;
template <typename stepper_t>
Acts::Propagator<stepper_t, Acts::Experimental::DetectorNavigator> makePropagator(
        std::shared_ptr<const Acts::Experimental::Detector> detector, std::shared_ptr<LUXEMagneticField::BField_t> BFieldPtr) {
    Acts::Experimental::DetectorNavigator::Config cfg{&(*detector)};
    cfg.resolvePassive = false;
    cfg.resolveMaterial = true;
    cfg.resolveSensitive = true;
    Acts::Experimental::DetectorNavigator navigator(
            cfg, Acts::getDefaultLogger("Detector Navigation", Acts::Logging::VERBOSE));

    stepper_t stepper(std::move(BFieldPtr));

    std::cout<<"Created stepper"<<std::endl;

    return Acts::Propagator<decltype(stepper), Acts::Experimental::DetectorNavigator>(
            std::move(stepper), std::move(navigator));
}

Acts::CurvilinearTrackParameters makeParameters(Acts::ActsScalar E) {
    Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity();
    // define a track in the transverse plane along x
    Acts::Vector4 mPos4(0., 0. , -0.618, 1_ns);
    return Acts::CurvilinearTrackParameters(mPos4, 90_degree, 90_degree,
                                            1_e / (E*1_GeV), cov, Acts::ParticleHypothesis::electron());
}

} // namespace LUXENavigator
