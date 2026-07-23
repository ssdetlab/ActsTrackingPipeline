#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"

#include <memory>

#include "TrackingPipeline/Simulation/IDigitizer.hpp"

/// @brief class constructing source links out of
/// simulated track parameters
class SimpleSourceLinkCreator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Digitizer
    std::shared_ptr<IDigitizer> hitDigitizer;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit SimpleSourceLinkCreator(const Config& cfg) : m_cfg(cfg) {}

  /// @brief generate Acts::SourceLink from the
  /// simulated track parameters
  ///
  /// @param gctx geometry context
  /// @param eventNumber current event number
  /// @param idx user-specified index
  /// @param trackParameters simulated track parameters
  /// @param rng random engine for sampling
  ///
  /// @return source link with the measurement
  Acts::SourceLink operator()(const Acts::GeometryContext& gctx,
                              std::size_t eventNumber, std::size_t idx,
                              const Acts::BoundTrackParameters& trackParameters,
                              RandomEngine& rng) const;

 private:
  /// Configuration struct
  Config m_cfg;
};
