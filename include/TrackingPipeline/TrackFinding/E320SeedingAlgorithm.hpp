#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/Infrastructure/DataHandle.hpp"
#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"
#include "TrackingPipeline/TrackFitting/FastGX2Fiter.hpp"

struct IdxTree {
  struct Node {
    Node() = delete;
    explicit Node(const std::pair<int, int>& data)
        : m_idx(data.first), m_geoId(data.second) {};

    int m_idx;
    int m_geoId;
    std::vector<std::shared_ptr<Node>> children;
  };

  using IdxContainer = std::vector<std::pair<int, int>>;

  IdxTree() = delete;

  IdxTree(const IdxContainer& container, const IdxContainer::iterator& root,
          const IdxContainer::iterator& rootEnd) {
    m_root = std::make_shared<Node>(*root);
    std::vector<std::shared_ptr<Node>> currentLayerNodes = {m_root};

    int layerId = rootEnd->second;
    IdxContainer layerIdxs;
    layerIdxs.reserve(container.size() / 5);
    for (auto it = rootEnd; it != container.end(); ++it) {
      int id = it->second;
      if (id == layerId) {
        layerIdxs.push_back(*it);
      } else {
        layerId = id;

        auto children = initNodes(layerIdxs);
        for (auto& node : currentLayerNodes) {
          addChildren(node, children);
        }
        currentLayerNodes = std::move(children);
        layerIdxs.clear();
        layerIdxs.push_back(*it);
        layerIdxs.reserve(container.size() / 5);
      }
      if (it == container.end() - 1) {
        auto children = initNodes(layerIdxs);
        for (auto& node : currentLayerNodes) {
          addChildren(node, children);
        }
      }
    }
  }

  std::vector<std::shared_ptr<Node>> initNodes(const IdxContainer& idxs) const {
    std::vector<std::shared_ptr<Node>> nodes;
    for (const auto& sl : idxs) {
      nodes.push_back(std::make_shared<Node>(sl));
    }
    return nodes;
  }

  std::vector<std::shared_ptr<Node>> addChildren(
      std::shared_ptr<Node>& parent,
      const std::vector<std::shared_ptr<Node>>& children) const {
    parent->children = children;
    return parent->children;
  }

  std::shared_ptr<Node> m_root;
};

class E320SeedingAlgorithm : public IAlgorithm {
 public:
  enum struct PropagationDirection : int { forward = 0, backward = 1 };

  /// @brief The nested configuration struct
  struct Config {
    /// HT seeder
    std::shared_ptr<HoughTransformSeeder> htSeeder;
    /// HT seeder options
    HoughTransformSeeder::Options htOptions;
    /// Input detector source links
    std::string inputSourceLinks;
    /// Input detector source links indices
    std::string inputDetSourceLinkIndices;
    /// Input BPM source links indices
    std::string inputBpmSourceLinkIndices;
    /// Output seeds
    std::string outputSeeds;
    /// Output track parameters
    std::string outputTrackParameters;
    /// Fast GX2 fitter
    std::shared_ptr<FastGX2Fitter> gx2Fitter;
    /// Number of GX2 fit iterations
    std::size_t nGX2Iterations;
    /// Maximum allowed chi2 of the GX2 fit
    double maxChi2;
    /// Reference surface
    const Acts::Surface* referenceSurface;
    /// Initial track state covariance prior
    Acts::BoundMatrix originCov;
    /// Lower cutoff on the number of layers in a seed
    std::size_t minLayers;
    /// Higher cutoff on the number of layers in a seed
    std::size_t maxLayers;
    /// Forward or backward propagation parameters
    PropagationDirection propDirection;
  };

  /// @brief Constructor
  E320SeedingAlgorithm(const Config& config, Acts::Logging::Level level);

  /// @brief Execute method
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  Acts::BoundMatrix transportCovToReference(
      const Acts::GeometryContext& gctx, const Acts::Vector3& refSurfacePoint,
      const Acts::Vector3& point, const Acts::Vector3& dir,
      const Acts::ActsSquareMatrix<6>& cov) const;

  double m_dipoleLength;
  double m_dipoleFieldStrength;

  /// Read data handles
  ReadDataHandle<std::vector<Acts::SourceLink>> m_inputSourceLinks{
      this, "InputSourceLinks"};

  ReadDataHandle<std::vector<std::size_t>> m_inputDetSourceLinkIndices{
      this, "InputDetSourceLinkIndices"};

  ReadDataHandle<std::vector<std::size_t>> m_inputBpmSourceLinkIndices{
      this, "InputBpmSourceLinkIndices"};

  /// Write data handles
  WriteDataHandle<IndexSeeds> m_outputSeeds{this, "OutputSeeds"};

  WriteDataHandle<std::vector<Acts::CurvilinearTrackParameters>>
      m_outputTrackParameters{this, "OutputTrackParameters"};
};
