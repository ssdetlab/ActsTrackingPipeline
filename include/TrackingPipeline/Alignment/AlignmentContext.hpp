#pragma once

#include <memory>
#include <stdexcept>

#include "TrackingPipeline/Alignment/AlignmentStore.hpp"

class AlignmentContext {
 public:
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
