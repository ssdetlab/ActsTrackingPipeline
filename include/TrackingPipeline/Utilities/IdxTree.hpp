#pragma once

#include <memory>
#include <utility>
#include <vector>

namespace detail {

class IdxTree {
 public:
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

  void constructTracks(const std::shared_ptr<IdxTree::Node>& root,
                       std::vector<std::size_t>& track,
                       std::vector<std::vector<std::size_t>>& tracks) const {
    track.push_back(root->m_idx);
    if (root->children.size() == 0) {
      tracks.push_back(track);
      track.pop_back();
      return;
    }
    for (auto& child : root->children) {
      constructTracks(child, track, tracks);
    }
    track.pop_back();
  }

  const std::shared_ptr<IdxTree::Node>& getRootNode() const { return m_root; }

 private:
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

}  // namespace detail
