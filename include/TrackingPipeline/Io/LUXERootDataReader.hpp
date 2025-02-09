#pragma once

#include "Acts/Definitions/Units.hpp"

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/LUXEGeometryConstraints.hpp"
#include "TrackingPipeline/Io/RootSimDataReader.hpp"

namespace LUXEROOTReader {

using namespace Acts::UnitLiterals;

/// @brief Global to local conversion
/// for the LUXE geometry
Acts::Vector2 convertToLoc(const Acts::Vector3& glob,
                           const Acts::GeometryIdentifier geoId,
                           const LUXEGeometry::GeometryOptions& gOpt) {
  // TODO: geoIds for the electron arm
  int nStave = geoId.sensitive() / 10 - 1;
  int nChip = geoId.sensitive() % 10 - 1;
  if (nStave % 2 == 0) {
    Acts::Vector2 loc = Acts::Vector2((glob.x() - gOpt.chipXEven.at(nChip)),
                                      (glob.y() - gOpt.chipY));
    return loc;
  } else {
    Acts::Vector2 loc = Acts::Vector2((glob.x() - gOpt.chipXOdd.at(nChip)),
                                      (glob.y() - gOpt.chipY));
    return loc;
  }
}

// Convert the momentum to the IP
// because Arka didn't make my life easy
TLorentzVector convertToIP(TLorentzVector mom,
                           double mass = 0.511 * Acts::UnitConstants::MeV) {
  TLorentzVector ipMom;
  ipMom.SetE(mom.E());
  ipMom.SetPx(0);
  ipMom.SetPy(mom.Py());
  ipMom.SetPz(
      std::sqrt(mom.E() * mom.E() - ipMom.Py() * ipMom.Py() - mass * mass));
  return ipMom;
}

/// @brief The ROOT file reader for the LUXE simulation
/// that knows about the true hits and the true momenta
///
/// @note Covariance is implemented as a diagonal matrix
/// of ALPIDE intrinsic resolutions
class LUXEROOTSimDataReader : public RootSimDataReader {};

auto defaultSimConfig() {
  LUXEROOTSimDataReader::Config config;
  return config;
}

}  // namespace LUXEROOTReader
