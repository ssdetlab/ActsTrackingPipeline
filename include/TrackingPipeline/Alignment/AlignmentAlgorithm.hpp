#pragma once

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>

#include <functional>
#include <memory>
#include <vector>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"

class AlignmentGroup {
 public:
  AlignmentGroup(const std::string& name,
                 const std::vector<Acts::GeometryIdentifier>& geoIds)
      : m_name(name), m_map(constructHierarchyMap(geoIds)) {}

  // Access the name of the group
  std::string getNameOfGroup() const { return m_name; }

  // Useful for testing
  bool has(Acts::GeometryIdentifier geoId) {
    auto it = m_map.find(geoId);
    return (it == m_map.end()) ? false : *it;
  }

 private:
  //  storing the name in the class
  std::string m_name;
  Acts::GeometryHierarchyMap<bool> m_map;

  Acts::GeometryHierarchyMap<bool> constructHierarchyMap(
      const std::vector<Acts::GeometryIdentifier>& geoIds) {
    std::vector<Acts::GeometryHierarchyMap<bool>::InputElement> ies;
    for (const auto& geoId : geoIds) {
      ies.emplace_back(geoId, true);
    }
    return Acts::GeometryHierarchyMap<bool>(ies);
  }
};

class AlignmentTransformUpdater {
 public:
  bool update(Acts::DetectorElementBase* detElement,
              const Acts::GeometryContext& gctx,
              const Acts::Transform3& transform) {
    return true;
  };
};

class AlignmentAlgorithm final : public IAlgorithm {
 public:
  using AlignmentResult = Acts::Result<ActsAlignment::AlignmentResult>;
  using AlignmentParameters =
      std::unordered_map<Acts::DetectorElementBase*, Acts::Transform3>;
  using TrackFitterOptions =
      Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory>;

  /// Alignment function that takes the above parameters and runs alignment
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class AlignmentFunction {
   public:
    virtual ~AlignmentFunction() = default;
    virtual AlignmentResult operator()(
        const std::vector<std::vector<Acts::SourceLink>>&,
        const std::vector<Acts::CurvilinearTrackParameters>&,
        const ActsAlignment::AlignmentOptions<TrackFitterOptions>&) = 0;
  };

  /// Create the alignment function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyway.
  static std::shared_ptr<AlignmentFunction> makeAlignmentFunction(
      std::shared_ptr<const Acts::Experimental::Detector> detector,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField);

  struct Config {
    /// Input track candidates
    std::string inputTrackCandidates;
    /// Output aligned parameters collection.
    std::string outputAlignmentParameters;
    /// Reference suface
    std::shared_ptr<Acts::Surface> referenceSurface;
    /// Type erased fitter function.
    std::shared_ptr<AlignmentFunction> align;
    /// The aligned transform updater
    ActsAlignment::AlignedTransformUpdater alignedTransformUpdater;
    /// The surfaces (with detector elements) to be aligned
    std::vector<Acts::DetectorElementBase*> alignedDetElements;
    /// KF options
    Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions;
    /// The alignment mask at each iteration
    std::map<unsigned int, std::bitset<6>> iterationState;
    /// Cutoff value for average chi2/ndf
    double chi2ONdfCutOff = 0.10;
    /// Cutoff value for delta of average chi2/ndf within a couple of iterations
    std::pair<std::size_t, double> deltaChi2ONdfCutOff = {10, 0.00001};
    /// Maximum number of iterations
    std::size_t maxNumIterations = 100;
    /// Number of tracks to be used for alignment
    int maxNumTracks = -1;
    std::vector<AlignmentGroup> m_groups;
  };

  /// Constructor of the alignment algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param level is the logging level
  AlignmentAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Framework execute method of the alignment algorithm
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ProcessCode execute(const AlgorithmContext& ctx) const override;

 private:
  Config m_cfg;

  ReadDataHandle<Seeds> m_inputTrackCandidates{this, "InputTrackCandidates"};

  WriteDataHandle<AlignmentParameters> m_outputAlignmentParameters{
      this, "OutputAlignmentParameters"};
};
