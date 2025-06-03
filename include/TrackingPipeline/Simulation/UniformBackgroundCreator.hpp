#pragma once

#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <utility>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Simulation/IMeasurementGenerator.hpp"

class UniformBackgroundCreator : public IMeasurementGenerator {
 public:
  using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                      Acts::Experimental::DetectorNavigator>;

  using TrackParameters = Acts::CurvilinearTrackParameters;

  /// @brief Nested configuration struct
  struct Config {
    std::size_t nMeasurements;
    std::vector<const Acts::Surface*> surfaces;
    std::pair<double, double> resolution;
  };

  /// @brief Constructor
  UniformBackgroundCreator(const Config& config) : m_cfg(config) {
    Acts::Vector2 stdDev = {m_cfg.resolution.first, m_cfg.resolution.second};
    m_cov = stdDev.cwiseProduct(stdDev).asDiagonal();
  };

  /// @brief Propagate track and create measurements
  std::tuple<std::vector<Acts::SourceLink>, SimClusters> gen(
      const AlgorithmContext& ctx, RandomEngine& rng,
      std::size_t id) const override {
    std::uniform_real_distribution uniform(0.0, 1.0);

    std::vector<Acts::SourceLink> sourceLinks;
    sourceLinks.reserve(m_cfg.surfaces.size() * m_cfg.nMeasurements);

    SimClusters clusters;
    clusters.reserve(m_cfg.surfaces.size() * m_cfg.nMeasurements);
    for (const auto* surf : m_cfg.surfaces) {
      const auto& bounds = surf->bounds().values();

      for (std::size_t i = 0; i < m_cfg.nMeasurements; i++) {
        double x = bounds.at(0) + uniform(rng) * (bounds.at(2) - bounds.at(0));
        double y = bounds.at(1) + uniform(rng) * (bounds.at(3) - bounds.at(1));

        SimpleSourceLink ssl(Acts::Vector2(x, y), m_cov, surf->geometryId(),
                             ctx.eventNumber, i);
        sourceLinks.push_back(Acts::SourceLink(ssl));

        SimCluster cluster{ssl, {}, false};
        clusters.push_back(cluster);
      }
    }
    return {sourceLinks, clusters};
  };

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  /// Covariance
  Acts::SquareMatrix2 m_cov;
};
