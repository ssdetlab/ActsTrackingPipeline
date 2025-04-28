#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include <Acts/Geometry/GeometryIdentifier.hpp>

#include <map>
#include <memory>

struct AlignmentContext {
  std::shared_ptr<std::map<Acts::GeometryIdentifier, Acts::Transform3>>
      alignmentStore = nullptr;

  /// Default constructor
  AlignmentContext() = default;

  /// Constructor with Store and context index
  AlignmentContext(
      std::shared_ptr<std::map<Acts::GeometryIdentifier, Acts::Transform3>>
          aStore)
      : alignmentStore(std::move(aStore)) {}
};
