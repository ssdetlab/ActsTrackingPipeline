#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Result.hpp"

class ITrackParametersEstimator {
 public:
  using Result = Acts::Result<Acts::CurvilinearTrackParameters>;

  virtual Result estimateParameters(
      const Acts::GeometryContext& gctx,
      const std::vector<Acts::SourceLink>& sourceLinks,
      const std::vector<std::size_t>& sourceLinkIndices,
      const Acts::Vector3& dir, const Acts::Vector3& point) const = 0;
};
