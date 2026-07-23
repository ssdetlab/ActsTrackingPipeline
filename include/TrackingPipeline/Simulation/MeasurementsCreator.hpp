#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Simulation/IMeasurementGenerator.hpp"
#include "TrackingPipeline/Simulation/IMomentumGenerator.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

using namespace Acts::UnitLiterals;

/// @brief Class creating measurements on the
/// sensitive surfaces crossed by a simulated track
class MeasurementsCreator : public IMeasurementGenerator {
 public:
  /// Propagator and track parameters shorthands
  using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                      Acts::Experimental::DetectorNavigator>;
  using TrackParameters = Acts::CurvilinearTrackParameters;

  /// Callable creating the source links out of true track parameters
  using SourceLinkCreator = Acts::Delegate<Acts::SourceLink(
      const Acts::GeometryContext& gctx, std::size_t eventNumber,
      std::size_t idx, const Acts::BoundTrackParameters& trackParameters,
      RandomEngine& rng)>;

  /// @brief Nested constraints struct
  struct Constraints {
    /// Surface frame of reference spatial cuts
    /// abotrting the propagation
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
    /// Map of surface ID to source link creator
    std::unordered_map<std::size_t, SourceLinkCreator> sourceLinkCreators;
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
  ///
  /// @param propagator propagator to be used
  /// @param cfg configuration struct
  MeasurementsCreator(const Propagator& propagator, const Config& cfg);

  /// @brief Propagate track and create measurements
  ///
  /// @param ctx current algorithm context
  /// @param rng random engine for sampling
  /// @param id user-specified source id
  /// @param sourceLinks container to store created source links
  /// @param simClusters container to store created sim clusters
  /// @param sourceLinksIndices container to store created source link indices
  ///
  /// @return number of created measurements
  std::size_t gen(const AlgorithmContext& ctx, RandomEngine& rng,
                  std::size_t id, std::vector<Acts::SourceLink>& sourceLinks,
                  SimClusters& simClusters,
                  std::vector<std::size_t>& sourceLinksIndices) const override;

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
