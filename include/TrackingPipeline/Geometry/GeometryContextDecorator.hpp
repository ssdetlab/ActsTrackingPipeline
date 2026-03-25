#pragma once

#include "Acts/Geometry/GeometryContext.hpp"

#include <string>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Alignment/AlignmentStore.hpp"
#include "TrackingPipeline/Infrastructure/IContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

class GeometryContextDecorator : public IContextDecorator {
 public:
  GeometryContextDecorator(
      const std::shared_ptr<AlignmentStore>& alignmentStorePtr)
      : m_alignmentStorePtr(alignmentStorePtr) {}

  ProcessCode decorate(AlgorithmContext& context) override {
    context.geoContext =
        Acts::GeometryContext{AlignmentContext(m_alignmentStorePtr)};
    return ProcessCode::SUCCESS;
  }

  const std::string& name() const override { return m_name; };

 private:
  std::shared_ptr<AlignmentStore> m_alignmentStorePtr = nullptr;

  std::string m_name = "GeometryContextDecorator";
};
