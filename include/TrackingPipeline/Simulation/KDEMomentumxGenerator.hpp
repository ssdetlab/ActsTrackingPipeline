#pragma once

#include <Acts/Definitions/Algebra.hpp>

#include <memory>
#include <vector>

#include "TrackingPipeline/Io/ITrackParamsReader.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/detail/NormalKDE.hpp"

/// @brief Class generating momentum vectors based on
/// kernel-density-estimation technique with user-provided
/// KDE sample
///
/// TODO: implement mean, covariance matrix, split away cpp
class KDEMomentumGenerator : public IMomentumGenerator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    std::shared_ptr<ITrackParamsReader> trackParamsReader;
    std::size_t nIterations;
    double sensitivity;
    Acts::Transform3 transform;
  };

  /// @brief Constructor
  ///
  /// @param config configuration struct
  explicit KDEMomentumGenerator(const Config& config) : m_cfg(config) {
    auto trackParams = m_cfg.trackParamsReader->read();

    std::vector<Acts::Vector3> sample;
    sample.reserve(trackParams.size());
    for (const auto& param : trackParams) {
      sample.emplace_back(
          Acts::Vector3{param.phi(), param.theta(), param.absoluteMomentum()});
    }
    m_kde = std::make_unique<NormalKDE<3>>(std::move(sample), m_cfg.nIterations,
                                           m_cfg.sensitivity);
  }

  /// @brief Generate momentum vector
  ///
  /// @param rng random engine for sampling
  ///
  /// @return Generated momentum vector
  Acts::Vector3 genMomentum(RandomEngine& rng) const override {
    Acts::Vector3 phiThetaE = m_kde->sample(rng);
    double phi = phiThetaE.x();
    double theta = phiThetaE.y();
    double E = phiThetaE.z();

    Acts::Vector3 dir{std::sin(theta) * std::cos(phi),
                      std::sin(theta) * std::sin(phi), std::cos(theta)};

    dir = m_cfg.transform * dir;

    return E * dir;
  }

 private:
  /// Configuration
  Config m_cfg;

  /// Normal KDE instance
  std::unique_ptr<NormalKDE<3>> m_kde;
};
