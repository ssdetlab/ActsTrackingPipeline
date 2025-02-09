#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Io/RootSimDataReader.hpp"

namespace E320Io {

using namespace Acts::UnitLiterals;

/// @brief The ROOT file reader for the LUXE simulation
/// that knows about the true trueHits and the true momenta
///
/// @note Covariance is implemented as a diagonal matrix
/// of ALPIDE intrinsic resolutions
class E320RootSimDataReader : public RootSimDataReader {
 public:
  /// @brief The configuration struct
  struct Config : public RootSimDataReader::Config {
    /// The geometry options
    E320Geometry::GeometryOptions gOpt;
  };

  E320RootSimDataReader(const Config& config, Acts::Logging::Level level)
      : RootSimDataReader(config, level), m_cfg(config) {
    m_actsToWorld = m_cfg.gOpt.actsToWorld.rotation().inverse();
  }

  std::string name() const override { return "E320RootSimDataReader"; }

 private:
  Config m_cfg;

  Acts::RotationMatrix3 m_actsToWorld;

  double m_pixSizeX = 27_um;
  double m_pixSizeY = 29_um;

  // Prepare the measurements
  // for the Sequencer pipeline
  void prepareMeasurements(const AlgorithmContext& context,
                           std::vector<Acts::SourceLink>* sourceLinks,
                           SimClusters* clusters) const override;
};

inline auto defaultSimConfig() {
  E320RootSimDataReader::Config config;
  config.treeName = "clusters";
  config.vVector3Keys = {"tru_hit", "tru_vertex"};
  config.vector3Keys = {"rglobal_geo"};
  config.vLorentzKeys = {"tru_p", "tru_p_ip"};
  config.vIntKeys = {"tru_trackId", "tru_parenttrackId", "tru_runId"};
  config.intKeys = {"eventId", "geoId", "xsize", "ysize", "size", "isSignal"};
  return config;
}

}  // namespace E320Io
