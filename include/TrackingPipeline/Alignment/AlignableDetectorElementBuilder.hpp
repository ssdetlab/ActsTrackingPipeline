#pragma once

#include "Acts/Detector/interface/IDetectorElementBuilder.hpp"

#include <memory>

#include "TrackingPipeline/Alignment/AlignableDetectorElement.hpp"

class AlignableDetectorElementBuilder
    : public Acts::Experimental::IDetectorElementBuilder {
 public:
  AlignableDetectorElementBuilder() = default;

  /// The interface method to be implemented by all detector
  /// component builder
  ///
  /// @param gctx The geometry context for this call
  ///
  /// @return an outgoing detector component
  std::vector<std::shared_ptr<Acts::DetectorElementBase>> construct(
      std::vector<std::shared_ptr<Acts::Surface>> surfaces,
      const Acts::GeometryContext& gctx) const override {
    std::vector<std::shared_ptr<Acts::DetectorElementBase>> detectorElements;
    detectorElements.reserve(surfaces.size());
    for (auto& surf : surfaces) {
      auto detElement = std::make_shared<AlignableDetectorElement>(
          surf, surf->transform(gctx));
      surf->assignDetectorElement(*detElement);
      detectorElements.push_back(detElement);
    }
    return detectorElements;
  }
};
