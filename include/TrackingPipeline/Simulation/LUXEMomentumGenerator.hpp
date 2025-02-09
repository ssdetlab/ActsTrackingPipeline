#pragma once

#include "Acts/Definitions/Units.hpp"

#include "TrackingPipeline/Geometry/LUXEGeometryConstraints.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"

namespace LUXESim {

using namespace Acts::UnitLiterals;

struct LUXEMomentumGenerator : public IMomentumGenerator {
  LUXEGeometry::GeometryOptions gOpt;

  double meanP = 2_MeV;
  double sigmaP = 1.8_MeV;

  double meanM = -2_MeV;
  double sigmaM = 1.8_MeV;

  double shapeZ = 3_GeV;
  double scaleZ = 1.2_GeV;

  Acts::Vector3 genMomentum(RandomEngine& rng) const override {
    std::normal_distribution<double> pDisP(meanP, sigmaP);
    std::normal_distribution<double> pDisM(meanM, sigmaM);
    std::gamma_distribution<double> pzDis(shapeZ, scaleZ);

    double px = (pDisP(rng) + pDisM(rng)) / 2;
    double py = (pDisP(rng) + pDisM(rng)) / 2;
    double pz = pzDis(rng) + 1;

    Acts::Vector3 momentum{px, py, pz};

    auto worldMomentum = gOpt.actsToWorld.rotation().inverse() * momentum;
    return worldMomentum;
  }
};

}  // namespace LUXESim
