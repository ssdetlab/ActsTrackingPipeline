#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>
#include <span>
#include <unordered_map>
#include <vector>

struct IdxTree {
  struct Node {
    Node() = delete;
    Node(std::pair<int, int> data) : m_idx(data.first), m_geoId(data.second) {};

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

struct TupleHash {
  std::size_t operator()(
      const std::tuple<std::uint16_t, std::uint16_t, std::uint16_t,
                       std::uint16_t>& t) const noexcept {
    return std::hash<long long>()(
        ((long long)std::get<0>(t) << 48) ^ ((long long)std::get<1>(t) << 32) ^
        ((long long)std::get<2>(t) << 16) ^ (long long)std::get<3>(t));
  }
};

class HoughTransformSeeder {
 public:
  using Cell =
      std::tuple<std::uint16_t, std::uint16_t, std::uint16_t, std::uint16_t>;
  using VotingMap = std::unordered_map<Cell, std::uint8_t, TupleHash>;

  using SourceLinkRef = std::reference_wrapper<const Acts::SourceLink>;

  struct Config {
    double boundBoxHalfX;
    double boundBoxHalfY;
    double boundBoxHalfZ;

    std::size_t nCellsThetaX;
    std::size_t nCellsRhoX;

    std::size_t nCellsThetaY;
    std::size_t nCellsRhoY;

    std::size_t minSeedSize;
    std::size_t maxSeedSize;

    std::size_t nLSIterations;
  };

  struct Options {
    double boundBoxCenterX;
    double boundBoxCenterY;
    double boundBoxCenterZ;

    int firstLayerId;
    int lastLayerId;
    int nLayers;

    int minCount;

    double maxChi2;
  };

  struct HTSeed {
    Acts::Vector3 lineRefPoint;
    Acts::Vector3 lineDir;
    std::vector<Acts::SourceLink> sourceLinks;
  };

  HoughTransformSeeder(const Config& cfg);

  std::vector<HTSeed> findSeeds(std::span<SourceLinkRef> sourceLinks,
                                const Options& opt) const;

 private:
  Config m_cfg;

  void fillVotingMap(VotingMap& votingMap, std::span<SourceLinkRef> points,
                     int nThreads, const Options& opt,
                     const Acts::Vector3& shift) const;

  std::vector<std::pair<int, int>> findLineSourceLinks(
      const std::span<SourceLinkRef>& sourceLinks,
      const Acts::Vector3& pointBL, const Acts::Vector3& dirBL,
      const Acts::Vector3& pointTR, const Acts::Vector3& dirTR,
      const Acts::Vector3& shift) const;

  std::vector<std::pair<int, int>> findLineSourceLinks(
      const std::span<SourceLinkRef>& sourceLinks, const Acts::Vector3& point,
      const Acts::Vector3& dir, double dist, const Acts::Vector3& shift) const;

  double orthogonalLeastSquares(const std::vector<SourceLinkRef>& sourceLinks,
                                Acts::Vector3& a, Acts::Vector3& b,
                                const Acts::Vector3& shift) const;

  double m_deltaThetaZ;
  double m_deltaRhoZ;

  double m_deltaThetaY;
  double m_deltaRhoY;

  double m_maxRhoY;
  double m_maxRhoZ;
};
