#pragma once

#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitter.hpp"
#include <Acts/Navigation/DetectorNavigator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/TrackFitting/GlobalChiSquareFitter.hpp>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"

class Gx2TrackFittingAlgorithm : public IAlgorithm {
 public:
  using ActionList = Acts::ActionList<>;
  using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

  using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                      Acts::Experimental::DetectorNavigator>;

  using TrackContainer = Acts::VectorTrackContainer;
  using Trajectory = Acts::VectorMultiTrajectory;

  /// @brief The nested configuration struct
  struct Config {
    /// The input collection
    std::string inputTrackCandidates;
    /// The output collection
    std::string outputTracks;
    /// fitter
    const Acts::Experimental::Gx2Fitter<Propagator, Trajectory>& fitter;
    /// Options
    Acts::Experimental::Gx2FitterOptions<Trajectory> options;
  };

  /// @brief Constructor
  Gx2TrackFittingAlgorithm(Config config, Acts::Logging::Level level)
      : IAlgorithm("TrackFittingAlgorithm", level), m_cfg(std::move(config)) {
    m_inputTrackCandidates.initialize(m_cfg.inputTrackCandidates);
    m_outputTracks.initialize(m_cfg.outputTracks);
  }
  ~Gx2TrackFittingAlgorithm() = default;

  /// @brief The execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<Seeds> m_inputTrackCandidates{this, "inputTrackCandidates"};

  WriteDataHandle<Tracks> m_outputTracks{this, "OutputTracks"};
};
