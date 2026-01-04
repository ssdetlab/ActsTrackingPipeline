#pragma once

#include "Acts/Geometry/GeometryContext.hpp"

#include <string>

#include "TrackingPipeline/Alignment/AlignmentContext.hpp"
#include "TrackingPipeline/Infrastructure/IContextDecorator.hpp"
#include "TrackingPipeline/Infrastructure/ProcessCode.hpp"

class GeometryContextDecorator : public IContextDecorator {
 public:
  GeometryContextDecorator(
      const std::shared_ptr<AlignmentContext::AlignmentStore>& alignmentStore)
      : m_alignmentStore(alignmentStore) {}

  ProcessCode decorate(AlgorithmContext& context) override {
    context.geoContext =
        Acts::GeometryContext{AlignmentContext(m_alignmentStore)};
    return ProcessCode::SUCCESS;
  }

  const std::string& name() const override { return m_name; };

 private:
  std::shared_ptr<AlignmentContext::AlignmentStore> m_alignmentStore = nullptr;

  std::string m_name = "GeometryContextDecorator";
};
