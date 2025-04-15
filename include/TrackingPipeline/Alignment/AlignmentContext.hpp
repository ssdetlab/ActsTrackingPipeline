#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cstddef>
#include <memory>

struct AlignmentContext {
  /// We have 2 different transforms
  std::shared_ptr<const std::array<Acts::Transform3, 2>> alignmentStore =
      nullptr;

  /// Context index
  std::size_t alignmentIndex{0};

  /// Default constructor
  AlignmentContext() = default;

  /// Constructor with Store and context index
  AlignmentContext(
      std::shared_ptr<const std::array<Acts::Transform3, 2>> aStore,
      unsigned int aIndex = 0)
      : alignmentStore(std::move(aStore)), alignmentIndex(aIndex) {}
};
