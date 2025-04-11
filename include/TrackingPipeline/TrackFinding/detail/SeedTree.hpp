#pragma once

#include "Acts/EventData/SourceLink.hpp"

#include <memory>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

struct SeedTree {
  struct Node {
    Node() = delete;
    Node(Acts::SourceLink sl) : m_sourceLink(sl) {};

    Acts::SourceLink m_sourceLink;
    std::vector<std::shared_ptr<Node>> children;
  };

  SeedTree() = delete;

  SeedTree(Seed seed) {
    root = std::make_shared<Node>(seed.sourceLinks.at(0));
    std::vector<std::shared_ptr<Node>> currentLayerNodes = {root};

    auto getLayerID = [](Acts::GeometryIdentifier geoId) {
      std::uint64_t id = geoId.sensitive();
      id /= 10;
      Acts::GeometryIdentifier layerId;
      layerId.setSensitive(id);
      return layerId;
    };

    Acts::GeometryIdentifier layerId =
        getLayerID(seed.sourceLinks.at(1).get<SimpleSourceLink>().geometryId());

    std::vector<Acts::SourceLink> layerSourceLinks;
    for (auto it = seed.sourceLinks.begin() + 1; it != seed.sourceLinks.end();
         ++it) {
      Acts::GeometryIdentifier id =
          getLayerID(it->get<SimpleSourceLink>().geometryId());
      if (id == layerId) {
        layerSourceLinks.push_back(*it);
      } else {
        layerId = id;

        auto children = initNodes(layerSourceLinks);
        for (auto& node : currentLayerNodes) {
          addChildren(node, children);
        }
        currentLayerNodes = children;

        layerSourceLinks.clear();
        layerSourceLinks.push_back(*it);
      }
      if (it == seed.sourceLinks.end() - 1) {
        auto children = initNodes(layerSourceLinks);
        for (auto& node : currentLayerNodes) {
          addChildren(node, children);
        }
      }
    }
  }

  std::vector<std::shared_ptr<Node>> initNodes(
      std::vector<Acts::SourceLink> m_sourceLink) {
    std::vector<std::shared_ptr<Node>> nodes;
    for (auto sl : m_sourceLink) {
      nodes.push_back(std::make_shared<Node>(sl));
    }
    return nodes;
  }

  std::vector<std::shared_ptr<Node>> addChildren(
      std::shared_ptr<Node> parent,
      std::vector<std::shared_ptr<Node>> children) {
    parent->children = children;
    return parent->children;
  }

  std::shared_ptr<Node> root;
};
