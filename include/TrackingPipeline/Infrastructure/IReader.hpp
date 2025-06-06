#pragma once

#include <utility>

#include "TrackingPipeline/Infrastructure/AlgorithmContext.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"
#include "TrackingPipeline/Infrastructure/SequenceElement.hpp"

/// Event data reader interface.
///
/// Read data from disk and add it to the event store. The reader can have
/// internal state and implementations are responsible to handle concurrent
/// calls.
class IReader : public SequenceElement {
 public:
  /// Provide range of available events or [0, SIZE_MAX) if undefined.
  ///
  /// The upper limit is exclusive, i.e. [0,3) means events 0, 1, and 2.
  virtual std::pair<std::size_t, std::size_t> availableEvents() const = 0;

  /// Read data for the requested event and write it into the event store.
  ///
  /// As a result of the parallelization and/or skipping events, this method
  /// will most likely not be called in order. Implementations must use the
  /// event number provided to select the proper data to be read.
  virtual ProcessCode read(const AlgorithmContext& context) = 0;

  /// Internal execute method forwards to the read method as mutable
  /// @param context The algorithm context
  ProcessCode internalExecute(const AlgorithmContext& context) final {
    return read(context);
  }

  /// Fulfill the algorithm interface
  ProcessCode initialize() override { return ProcessCode::SUCCESS; }

  /// Fulfill the algorithm interface
  ProcessCode finalize() override { return ProcessCode::SUCCESS; }
};
