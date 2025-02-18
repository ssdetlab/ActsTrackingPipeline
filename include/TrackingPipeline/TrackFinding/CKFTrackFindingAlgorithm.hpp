#pragma once

#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/Navigation/DetectorNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include <Acts/Utilities/Holders.hpp>

#include <cstddef>
#include <memory>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"

/// The map(-like) container accessor
template <typename container_t>
struct SimpleContainerAccessor {
  using Container = container_t;
  using Key = typename container_t::key_type;
  using Value = typename container_t::mapped_type;

  /// This iterator adapter is needed to have the deref operator return a single
  /// source link instead of the map pair <GeometryIdentifier,SourceLink>
  struct Iterator {
    using BaseIterator = typename container_t::const_iterator;

    using iterator_category = typename BaseIterator::iterator_category;
    using value_type = typename BaseIterator::value_type;
    using difference_type = typename BaseIterator::difference_type;
    using pointer = typename BaseIterator::pointer;
    using reference = typename BaseIterator::reference;

    Iterator& operator++() {
      ++m_iterator;
      return *this;
    }

    bool operator==(const Iterator& other) const {
      return m_iterator == other.m_iterator;
    }

    bool operator!=(const Iterator& other) const { return !(*this == other); }

    Acts::SourceLink operator*() const {
      const auto& sl = m_iterator->second;
      return Acts::SourceLink{sl};
    }

    BaseIterator m_iterator;
  };

  // pointer to the container
  const Container* container = nullptr;

  // get the range of elements with requested key
  std::pair<Iterator, Iterator> range(const Acts::Surface& surface) const {
    assert(container != nullptr);
    auto [begin, end] = container->equal_range(surface.geometryId());
    return {Iterator{begin}, Iterator{end}};
  }
};

class CKFTrackFindingAlgorithm : public IAlgorithm {
 public:
  /// Propagator definitions
  using ActionList = Acts::ActionList<>;
  using AbortList = Acts::AbortList<Acts::EndOfWorldReached>;

  using Propagator = Acts::Propagator<Acts::EigenStepper<>,
                                      Acts::Experimental::DetectorNavigator>;
  /// Track containers
  using TrackContainer = Acts::TrackContainer<Acts::VectorTrackContainer,
                                              Acts::VectorMultiTrajectory,
                                              std::shared_ptr>;
  using TrackStateContainerBackend =
      typename TrackContainer::TrackStateContainerBackend;

  /// Sourcelink containers
  using SimpleSourceLinkContainer =
      std::unordered_multimap<Acts::GeometryIdentifier, SimpleSourceLink>;
  using SimpleSourceLinkAccessor =
      SimpleContainerAccessor<SimpleSourceLinkContainer>;

  /// Options
  using CKFOptions =
      Acts::CombinatorialKalmanFilterOptions<SimpleSourceLinkAccessor::Iterator,
                                             TrackContainer>;

  /// @brief The nested configuration struct
  struct Config {
    /// CKF finder
    const Acts::CombinatorialKalmanFilter<Propagator, TrackContainer>& ckf;
    /// CKF extensions
    Acts::CombinatorialKalmanFilterExtensions<TrackContainer> extensions;
    /// The input collection
    std::string inputSeeds;
    /// The output collection before
    std::string outputTrackCandidates;
    /// The output collection before
    std::string outputTrackView;
    /// Minimum number of source links
    std::size_t minCandidateSize;
    /// Maximum number of source links
    std::size_t maxCandidateSize;
    /// Maximum number of steps for the propagator
    std::size_t maxSteps;
  };

  /// @brief Constructor
  CKFTrackFindingAlgorithm(Config config, Acts::Logging::Level level)
      : IAlgorithm("CKFTrackFindingAlgorithm", level),
        m_cfg(std::move(config)) {
    m_inputSeeds.initialize(m_cfg.inputSeeds);
    m_outputTrackCandidates.initialize(m_cfg.outputTrackCandidates);
    m_outputTrackView.initialize(m_cfg.outputTrackView);
  }
  ~CKFTrackFindingAlgorithm() = default;

  /// @brief The execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override;

 private:
  Config m_cfg;

  ReadDataHandle<Seeds> m_inputSeeds{this, "InputSeeds"};

  WriteDataHandle<Seeds> m_outputTrackCandidates{this, "OutputTrackCandidates"};
  WriteDataHandle<Tracks> m_outputTrackView{this, "OutputTrackView"};
};
