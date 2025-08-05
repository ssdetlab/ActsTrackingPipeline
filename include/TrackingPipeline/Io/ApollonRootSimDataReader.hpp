#pragma once

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/Io/RootSimDataReader.hpp"

namespace ApollonIo {

using namespace Acts::UnitLiterals;

/// @brief The ROOT file reader for the Apollon simulation
/// that knows about the true trueHits and the true momenta
///
/// @note Covariance is implemented as a diagonal matrix
/// of ALPIDE intrinsic resolutions
class ApollonRootSimDataReader : public RootSimDataReader {
 public:
  /// @brief The configuration struct
  struct Config : public RootSimDataReader::Config {
    /// The geometry options
    ApollonGeometry::GeometryOptions gOpt;

    Acts::SourceLinkSurfaceAccessor surfaceAccessor;
  };

  ApollonRootSimDataReader(const Config& config, Acts::Logging::Level level)
      : RootSimDataReader(config, level), m_cfg(config) {}

  std::string name() const override { return "ApollonRootSimDataReader"; }

 private:
  Config m_cfg;

  // Prepare the measurements
  // for the Sequencer pipeline
  void prepareMeasurements(const AlgorithmContext& context,
                           std::vector<Acts::SourceLink>* sourceLinks,
                           SimClusters* clusters) const override;
};

inline auto defaultSimConfig() {
  ApollonRootSimDataReader::Config config;
  config.treeName = "particles";
  config.vVector3Keys = {"hitPosGlobal", "ipMomDir", "vertex", "hitMomDir"};
  config.vVector2Keys = {"hitPosLocal"};
  config.vIntKeys = {"geoId", "isSignal", "trackId", "parentTrackId", "runId"};
  config.intKeys = {"eventId"};
  config.vDoubleKeys = {"ipE", "hitE"};
  return config;
}

}  // namespace ApollonIo
