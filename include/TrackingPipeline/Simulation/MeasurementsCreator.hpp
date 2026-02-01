#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include <memory>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Simulation/IDigitizer.hpp"
#include "TrackingPipeline/Simulation/IMeasurementGenerator.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

using namespace Acts::UnitLiterals;

/// @brief Class creating measurements along
/// simulated tracks
///
/// Algorithm propagates tracks thorugh the geometry and
/// create measurements on the encountered senstive
/// surfaces
class MeasurementsCreator : public IMeasurementGenerator {
 public:
  using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                      Acts::Experimental::DetectorNavigator>;

  using TrackParameters = Acts::CurvilinearTrackParameters;

  struct Constraints {
    double minLocX;
    double maxLocX;

    double minLocY;
    double maxLocY;
  };

  /// @brief Nested configuration struct
  struct Config {
    /// Vertex generator
    std::shared_ptr<IVertexGenerator> vertexGenerator;
    /// Momentum generator
    std::shared_ptr<IMomentumGenerator> momentumGenerator;
    /// Digitizer
    std::shared_ptr<IDigitizer> hitDigitizer;
    /// Reference surface
    const Acts::Surface* referenceSurface;
    /// Maximum number of steps
    /// to propagate
    std::size_t maxSteps;
    /// Is signal flag
    bool isSignal;
    /// Particle hypothesis
    Acts::ParticleHypothesis hypothesis;
    /// Particle charge
    double charge;
    /// Constraint surfaces cuts
    std::unordered_map<Acts::GeometryIdentifier, Constraints> constraints;
  };

  /// @brief Constructor
  MeasurementsCreator(const Propagator propagator, const Config& config);

  /// @brief Propagate track and create measurements
  std::tuple<std::vector<Acts::SourceLink>, SimClusters> gen(
      const AlgorithmContext& ctx, RandomEngine& rng,
      std::size_t id) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  /// IP process covariance
  Acts::FreeMatrix m_freeIpCov;

  /// Propagator instance
  Propagator m_propagator;
};
