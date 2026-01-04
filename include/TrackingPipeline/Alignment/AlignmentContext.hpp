#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <memory>
#include <stdexcept>
#include <unordered_map>

class AlignmentContext {
 public:
  using AlignmentStore =
      std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3>;

  /// Default constructor
  AlignmentContext() = delete;

  /// Constructor with Store and context index
  AlignmentContext(const std::shared_ptr<AlignmentStore>& aStore) {
    if (aStore == nullptr) {
      throw std::runtime_error("Invalid alignment store initialization");
    }
    m_alignmentStore = aStore;
  }

  const AlignmentStore& alignmentStore() const { return *m_alignmentStore; }

  AlignmentStore& alignmentStore() { return *m_alignmentStore; }

  const std::shared_ptr<AlignmentStore> alignmentStorePtr() const {
    return m_alignmentStore;
  }

  std::shared_ptr<AlignmentStore> alignmentStorePtr() {
    return m_alignmentStore;
  }

 private:
  std::shared_ptr<AlignmentStore> m_alignmentStore = nullptr;
};
